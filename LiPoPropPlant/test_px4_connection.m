% TEST_PX4_CONNECTION  Standalone test – verify TCP + MAVLink with PX4 SITL
%
% Run this in MATLAB *before* touching Simulink.  It opens a TCP server
% on port 4560, waits for PX4 SITL to connect, exchanges a few frames,
% and prints what it receives.
%
% USAGE:
%   1. Run this script in MATLAB on Windows.
%   2. In WSL2 Ubuntu, launch PX4:
%        PX4_SIM_HOST_ADDR=<windows_ip> SYS_AUTOSTART=10016 \
%            PX4_SIM_SPEED_FACTOR=1 make px4_sitl_default none_iris
%   3. Watch the MATLAB command window for connection + motor data.
%
% Press Ctrl+C to stop.
% ---------------------------------------------------------------

fprintf('\n===== PX4 TCP/MAVLink Connection Test =====\n\n');

% ---- Create TCP server ------------------------------------------
try
    srv = java.net.ServerSocket();
    srv.setReuseAddress(true);
    srv.bind(java.net.InetSocketAddress(4560));
    srv.setSoTimeout(500);
    fprintf('Listening on port 4560 ... start PX4 SITL now.\n');
catch e
    fprintf('ERROR: cannot bind port 4560: %s\n', e.message);
    fprintf('  → Close other MATLAB instances or Simulink models first.\n');
    return;
end

% ---- Wait for PX4 to connect ------------------------------------
client = [];
while isempty(client)
    try
        client = srv.accept();
        fprintf('PX4 connected from %s\n', char(client.getRemoteSocketAddress()));
    catch
        fprintf('.');
    end
end

client.setTcpNoDelay(true);
client.setSoTimeout(5);
in_s  = client.getInputStream();
out_s = client.getOutputStream();

% ---- Helper: send bytes ------------------------------------------
send = @(bytes) send_bytes(out_s, uint8(bytes));

% ---- Main loop: send fake sensor data, read actuator commands ----
seq = uint8(0);
sim_time_us = uint64(0);
dt_us = uint64(4000);          % 250 Hz

fprintf('\nSending HIL_SENSOR + HIL_GPS ... reading motor commands.\n');
fprintf('(Ctrl+C to stop)\n\n');

rx_buf = uint8([]);
step = 0;

cleanupObj = onCleanup(@() cleanup_sockets(srv, client));

while true
    step = step + 1;
    sim_time_us = sim_time_us + dt_us;
    
    % ---- Send HEARTBEAT every 250 steps (~1 Hz) -----------------
    if mod(step, 250) == 1
        hb_payload = heartbeat_payload();
        [frame, seq] = mav2_frame(0, hb_payload, 50, seq);
        send(frame);
    end
    
    % ---- Send HIL_SENSOR every step (250 Hz) ---------------------
    % Fake data: sitting on ground, 1g down, Abu Dhabi mag field
    sensor_pl = hil_sensor_payload( ...
        sim_time_us, ...
        single(0), single(0), single(-9.80665), ...  % accel FRD
        single(0), single(0), single(0), ...          % gyro
        single(0.23), single(0.04), single(0.40), ... % mag (Gauss)
        single(1013.25), single(0), ...               % abs/diff press
        single(0), single(25), ...                    % alt, temp
        uint32(hex2dec('1FFF')));
    [frame, seq] = mav2_frame(107, sensor_pl, 108, seq);
    send(frame);
    
    % ---- Send HIL_GPS every 25 steps (~10 Hz) --------------------
    if mod(step, 25) == 0
        gps_pl = hil_gps_payload( ...
            sim_time_us, ...
            int32(24.4539 * 1e7), ...    % lat degE7 (Abu Dhabi)
            int32(54.3773 * 1e7), ...    % lon degE7
            int32(0), ...                % alt mm (sea level)
            uint16(100), uint16(100), ...% eph/epv cm
            uint16(0), ...               % ground speed cm/s
            int16(0), int16(0), int16(0), ... % vn/ve/vd cm/s
            uint16(0), ...               % COG cdeg
            uint8(3), uint8(12));        % 3D fix, 12 sats
        [frame, seq] = mav2_frame(113, gps_pl, 124, seq);
        send(frame);
    end
    
    % ---- Read from PX4 -------------------------------------------
    try
        avail = in_s.available();
        if avail > 0
            new_bytes = uint8(zeros(1, avail));
            for k = 1:avail
                b = in_s.read();
                if b < 0; break; end
                new_bytes(k) = uint8(b);
            end
            rx_buf = [rx_buf, new_bytes]; %#ok<AGROW>
        end
    catch
        % no data yet
    end
    
    % ---- Parse MAVLink frames ------------------------------------
    while length(rx_buf) >= 12
        idx = find(rx_buf == hex2dec('FD'), 1);
        if isempty(idx); rx_buf = uint8([]); break; end
        if idx > 1; rx_buf = rx_buf(idx:end); end
        if length(rx_buf) < 2; break; end
        
        plen = double(rx_buf(2));
        flen = 10 + plen + 2;
        if length(rx_buf) < flen; break; end
        
        msg_id = uint32(rx_buf(8)) + ...
                 bitshift(uint32(rx_buf(9)), 8) + ...
                 bitshift(uint32(rx_buf(10)), 16);
        payload = rx_buf(11 : 10+plen);
        rx_buf = rx_buf(flen+1:end);
        
        if msg_id == 93 && plen >= 72
            % HIL_ACTUATOR_CONTROLS
            t = zeros(1,4);
            for m = 1:4
                off = 8 + (m-1)*4;
                t(m) = double(typecast(uint8(payload(off+1:off+4)), 'single'));
            end
            t = max(0, min(1, t));
            fprintf('t=%.2fs  Motors: [%.4f  %.4f  %.4f  %.4f]\n', ...
                double(sim_time_us)/1e6, t(1), t(2), t(3), t(4));
        elseif msg_id == 0
            fprintf('t=%.2fs  Received HEARTBEAT from PX4\n', ...
                double(sim_time_us)/1e6);
        else
            % Uncomment to see all messages:
            % fprintf('t=%.2fs  msg_id=%d  len=%d\n', ...
            %     double(sim_time_us)/1e6, msg_id, plen);
        end
    end
    
    % Pace at roughly real-time (250 Hz → 4 ms per step)
    pause(0.004);
end

% =====================================================================
function cleanup_sockets(srv, cli)
    fprintf('\nCleaning up sockets...\n');
    try cli.close(); catch; end
    try srv.close(); catch; end
    fprintf('Done.\n');
end

function ok = send_bytes(out_s, data)
    ok = true;
    try
        signed = int8(data);
        signed(data > 127) = int8(double(data(data > 127)) - 256);
        out_s.write(signed);
        out_s.flush();
    catch
        ok = false;
    end
end

function [frame, seq] = mav2_frame(msg_id, payload, crc_extra, seq)
    SYS_ID = uint8(142); COMP_ID = uint8(142);
    len = uint8(length(payload));
    id0 = uint8(bitand(msg_id, 255));
    id1 = uint8(bitand(bitshift(msg_id, -8), 255));
    id2 = uint8(bitand(bitshift(msg_id, -16), 255));
    hdr = uint8([len, 0, 0, seq, SYS_ID, COMP_ID, id0, id1, id2]);
    crc = crc16([hdr, payload, crc_extra]);
    frame = uint8([hex2dec('FD'), hdr, payload, ...
                   uint8(bitand(crc,255)), uint8(bitand(bitshift(crc,-8),255))]);
    seq = uint8(mod(double(seq)+1, 256));
end

function crc = crc16(data)
    crc = uint16(hex2dec('FFFF'));
    for i = 1:length(data)
        tmp = bitxor(uint8(data(i)), uint8(bitand(crc, uint16(255))));
        tmp = bitxor(tmp, uint8(bitand(bitshift(uint16(tmp),4), uint16(255))));
        t16 = uint16(tmp);
        crc = bitxor(bitshift(crc,-8), bitxor(bitshift(t16,8), ...
              bitxor(bitshift(t16,3), bitshift(t16,-4))));
    end
end

function pl = heartbeat_payload()
    pl = uint8([typecast(uint32(0),'uint8'), 6, 8, 0, 4, 3]);
end

function pl = hil_sensor_payload(t, ax,ay,az, gx,gy,gz, mx,my,mz, ...
                                  ap,dp, pa, temp, fu)
    pl = uint8(zeros(1,64));
    pl(1:8)=typecast(t,'uint8');
    pl(9:12)=typecast(ax,'uint8');  pl(13:16)=typecast(ay,'uint8');
    pl(17:20)=typecast(az,'uint8'); pl(21:24)=typecast(gx,'uint8');
    pl(25:28)=typecast(gy,'uint8'); pl(29:32)=typecast(gz,'uint8');
    pl(33:36)=typecast(mx,'uint8'); pl(37:40)=typecast(my,'uint8');
    pl(41:44)=typecast(mz,'uint8'); pl(45:48)=typecast(ap,'uint8');
    pl(49:52)=typecast(dp,'uint8'); pl(53:56)=typecast(pa,'uint8');
    pl(57:60)=typecast(temp,'uint8'); pl(61:64)=typecast(fu,'uint8');
end

function pl = hil_gps_payload(t, lat,lon,alt, eph,epv, vel, vn,ve,vd, cog, fix, sat)
    pl = uint8(zeros(1,36));
    pl(1:8)=typecast(t,'uint8');
    pl(9:12)=typecast(lat,'uint8');   pl(13:16)=typecast(lon,'uint8');
    pl(17:20)=typecast(alt,'uint8');  pl(21:22)=typecast(eph,'uint8');
    pl(23:24)=typecast(epv,'uint8');  pl(25:26)=typecast(vel,'uint8');
    pl(27:28)=typecast(vn,'uint8');   pl(29:30)=typecast(ve,'uint8');
    pl(31:32)=typecast(vd,'uint8');   pl(33:34)=typecast(cog,'uint8');
    pl(35)=fix; pl(36)=sat;
end