function throttles = px4_sim_bridge_V1_old(sensor_data)
% PX4_SIM_BRIDGE  Simulink ←→ PX4 SITL bridge via MAVLink v2 over TCP
%
% Drop-in replacement for 'px4_tcp_server'. Use inside an Interpreted
% MATLAB Function block.  Mux ALL plant outputs into one vector and feed
% it here; the function handles MAVLink encoding, TCP transport, and
% MAVLink decoding internally.
%
% ── INPUT (sensor_data) ──────────────────────────────────────────────
%  Index   Quantity                  Units       Frame
%  1-3     ax, ay, az               m/s²        body (FRD)
%  4-6     gx, gy, gz               rad/s       body (FRD)
%  7-9     mx, my, mz               Gauss       body (FRD)
%  10      abs_pressure             mbar (hPa)
%  11      diff_pressure            mbar (hPa)   (set 0 if unused)
%  12      pressure_alt             m            barometric altitude
%  13      temperature              °C
%  14      lat                      deg          WGS-84
%  15      lon                      deg          WGS-84
%  16      alt_msl                  m            above MSL
%  17-19   vn, ve, vd               m/s          NED
%  20      sim_time                 s            simulation clock
%
% ── OUTPUT (throttles) ──────────────────────────────────────────────
%  [throttle_1, throttle_2, throttle_3, throttle_4]   normalised [0,1]
%
% ── QUICK WIRING IN SIMULINK ───────────────────────────────────────
%  1. Mux your 6-DOF + environment outputs into the 20-element vector
%     described above.
%  2. Connect the Mux to an "Interpreted MATLAB Function" block.
%     • MATLAB function:  px4_sim_bridge
%     • Output dimensions: 4
%     • Output signal type: double
%  3. Demux the output into your four motor throttle channels.
%  4. Solver: Fixed-step, 4 ms (250 Hz) – matches PX4 IMU rate.
%  5. Start the Simulink simulation FIRST, then launch PX4 SITL.
%
% ── NOTES ──────────────────────────────────────────────────────────
%  • PX4 SITL is the TCP CLIENT; this function is the TCP SERVER on
%    port 4560.  PX4 connects to us.
%  • Lockstep: PX4 advances its clock only when it receives
%    HIL_SENSOR with a new time_usec.  We drive time from sim_time.
%  • HIL_GPS is sent at ~10 Hz (every 25th step at 250 Hz).
%  • A HEARTBEAT is sent at ~1 Hz.
%  • Frame convention: PX4 uses NED / FRD.  Make sure your 6-DOF
%    block outputs accelerations and angular rates in body-FRD, and
%    velocity in NED.  If your model uses NWU/FLU, negate y and z
%    components before muxing.
%
% Pavel's custom PX4 v1.17.0 / Simulink 2022b HIL simulation
% ---------------------------------------------------------------

persistent server_sock client_sock in_strm out_strm ...
           is_connected seq_cnt gps_divider hb_divider ...
           rx_buf last_throttles init_flag

% ----- Output defaults -------------------------------------------
throttles = zeros(1, 4);

% ----- First-call initialisation ----------------------------------
if isempty(init_flag)
    server_sock   = [];
    client_sock   = [];
    in_strm       = [];
    out_strm      = [];
    is_connected  = false;
    seq_cnt       = uint8(0);
    gps_divider   = uint32(0);
    hb_divider    = uint32(0);
    rx_buf        = uint8([]);
    last_throttles = zeros(1,4);
    init_flag     = true;
end

% ----- Step 1: Create TCP server if not yet open -------------------
if isempty(server_sock)
    try
        server_sock = java.net.ServerSocket();
        server_sock.setReuseAddress(true);
        server_sock.bind(java.net.InetSocketAddress(4560));
        server_sock.setSoTimeout(1);           % 1 ms non-blocking
        fprintf('[px4_sim_bridge] TCP server listening on port 4560\n');
    catch e
        fprintf('[px4_sim_bridge] Port 4560 busy: %s\n', e.message);
        server_sock = [];
        return;
    end
end

% ----- Step 2: Accept PX4 connection (non-blocking) ---------------
if ~is_connected
    try
        client_sock = server_sock.accept();     % blocks ≤ 1 ms
        client_sock.setTcpNoDelay(true);
        client_sock.setSoTimeout(1);
        in_strm  = client_sock.getInputStream();
        out_strm = client_sock.getOutputStream();
        is_connected = true;
        fprintf('[px4_sim_bridge] *** PX4 connected! ***\n');
    catch
        % Timeout – PX4 not connected yet, return silently
        return;
    end
end

% Return last known throttles while we process this step
throttles = last_throttles;

% ----- Step 3: Unpack the input vector ----------------------------
ax   = sensor_data(1);   ay   = sensor_data(2);   az   = sensor_data(3);
gx   = sensor_data(4);   gy   = sensor_data(5);   gz   = sensor_data(6);
mx   = sensor_data(7);   my   = sensor_data(8);   mz   = sensor_data(9);
abs_press  = sensor_data(10);
diff_press = sensor_data(11);
press_alt  = sensor_data(12);
temp_c     = sensor_data(13);
lat_deg    = sensor_data(14);
lon_deg    = sensor_data(15);
alt_msl    = sensor_data(16);
vn = sensor_data(17);  ve = sensor_data(18);  vd = sensor_data(19);
sim_time_s = sensor_data(20);

time_usec = uint64(sim_time_s * 1e6);

% ----- Step 4: Send HEARTBEAT at ~1 Hz ----------------------------
hb_divider = hb_divider + 1;
if hb_divider >= 250                        % every 250 steps @ 250 Hz
    hb_divider = uint32(0);
    hb_payload = build_heartbeat_payload();
    [hb_frame, seq_cnt] = mavlink2_frame(0, hb_payload, 50, seq_cnt);
    tcp_send(out_strm, hb_frame);
end

% ----- Step 5: Send HIL_SENSOR every step (250 Hz) ----------------
fields_updated = uint32(hex2dec('1FFF'));    % all fields valid
sensor_payload = build_hil_sensor_payload( ...
    time_usec, ...
    single(ax), single(ay), single(az), ...
    single(gx), single(gy), single(gz), ...
    single(mx), single(my), single(mz), ...
    single(abs_press), single(diff_press), ...
    single(press_alt), single(temp_c), ...
    fields_updated);
[sensor_frame, seq_cnt] = mavlink2_frame(107, sensor_payload, 108, seq_cnt);
if ~tcp_send(out_strm, sensor_frame)
    % Connection lost
    handle_disconnect();
    is_connected = false;
    return;
end

% ----- Step 6: Send HIL_GPS at ~10 Hz -----------------------------
gps_divider = gps_divider + 1;
if gps_divider >= 25                        % every 25 steps @ 250 Hz
    gps_divider = uint32(0);
    
    % Ground speed
    v_ground = sqrt(vn^2 + ve^2);
    % Course over ground (deg, 0=North, CW positive)
    cog_deg = atan2d(ve, vn);
    if cog_deg < 0; cog_deg = cog_deg + 360; end
    
    gps_payload = build_hil_gps_payload( ...
        time_usec, ...
        int32(lat_deg * 1e7), ...       % degE7
        int32(lon_deg * 1e7), ...       % degE7
        int32(alt_msl * 1e3), ...       % mm
        uint16(100), ...                % eph = 1.00 m (cm)
        uint16(100), ...                % epv = 1.00 m (cm)
        uint16(v_ground * 100), ...     % cm/s
        int16(vn * 100), ...            % cm/s
        int16(ve * 100), ...
        int16(vd * 100), ...
        uint16(cog_deg * 100), ...      % cdeg
        uint8(3), ...                   % fix_type = 3D fix
        uint8(12));                     % satellites
    [gps_frame, seq_cnt] = mavlink2_frame(113, gps_payload, 124, seq_cnt);
    tcp_send(out_strm, gps_frame);
end

% ----- Step 7: Read data from PX4 ---------------------------------
try
    avail = in_strm.available();
    if avail > 0
        new_bytes = uint8(zeros(1, avail));
        for k = 1:avail
            b = in_strm.read();
            if b < 0; break; end        % EOF
            new_bytes(k) = uint8(b);
        end
        rx_buf = [rx_buf, new_bytes];
    end
catch
    handle_disconnect();
    is_connected = false;
    return;
end

% ----- Step 8: Parse MAVLink frames from rx buffer -----------------
while length(rx_buf) >= 12               % minimum MAVLink v2 frame
    % Find MAVLink v2 start byte
    idx = find(rx_buf == hex2dec('FD'), 1);
    if isempty(idx)
        rx_buf = uint8([]);
        break;
    end
    if idx > 1
        rx_buf = rx_buf(idx:end);         % discard leading junk
    end
    
    if length(rx_buf) < 2; break; end
    payload_len = rx_buf(2);
    frame_len   = 10 + double(payload_len) + 2;  % header(10) + payload + crc(2)
    
    if length(rx_buf) < frame_len
        break;                            % incomplete frame, wait
    end
    
    % Extract message ID (3 bytes, little-endian)
    msg_id = uint32(rx_buf(8)) + ...
             bitshift(uint32(rx_buf(9)), 8) + ...
             bitshift(uint32(rx_buf(10)), 16);
    
    % Payload bytes
    payload = rx_buf(11 : 10 + double(payload_len));
    
    % TODO: CRC check (skipped for speed; PX4 is trusted source)
    
    % Parse HIL_ACTUATOR_CONTROLS (msg ID 93)
    if msg_id == 93 && payload_len >= 72
        % controls[0..3] are at payload bytes 9..24 (after 8-byte time_usec)
        for m = 1:4
            offset = 8 + (m-1)*4;        % 0-indexed: 8,12,16,20
            throttles(m) = double(bytes_to_single( ...
                payload(offset+1 : offset+4)));
        end
        % Clamp to [0, 1]
        throttles = max(0, min(1, throttles));
        last_throttles = throttles;
    end
    
    % Consume this frame
    rx_buf = rx_buf(frame_len+1 : end);
end

% === Nested helper functions ======================================

    function handle_disconnect()
        fprintf('[px4_sim_bridge] Connection lost – waiting for reconnect\n');
        try client_sock.close(); catch; end
        client_sock = []; in_strm = []; out_strm = [];
    end

end  % px4_sim_bridge

% =====================================================================
%  LOCAL FUNCTIONS  (accessible only from this file)
% =====================================================================

function ok = tcp_send(out_strm, data)
% Send uint8 array over Java output stream.
    ok = true;
    try
        % Java expects signed bytes: convert uint8 → int8
        signed = int8(data);
        signed(data > 127) = int8(double(data(data > 127)) - 256);
        out_strm.write(signed);
        out_strm.flush();
    catch
        ok = false;
    end
end

function [frame, seq] = mavlink2_frame(msg_id, payload, crc_extra, seq)
% Build a complete MAVLink v2 frame.
%   msg_id    : uint32 message ID
%   payload   : uint8 row vector
%   crc_extra : uint8 message-specific CRC seed
%   seq       : uint8 sequence counter (updated on return)
    SYS_ID  = uint8(142);   % Simulator convention
    COMP_ID = uint8(142);
    
    len = uint8(length(payload));
    
    % Message ID as 3 little-endian bytes
    id0 = uint8(bitand(msg_id, 255));
    id1 = uint8(bitand(bitshift(msg_id, -8), 255));
    id2 = uint8(bitand(bitshift(msg_id, -16), 255));
    
    header = uint8([len, 0, 0, seq, SYS_ID, COMP_ID, id0, id1, id2]);
    
    % X.25 CRC over header + payload + crc_extra
    crc = mavlink_crc16([header, payload, crc_extra]);
    crc_lo = uint8(bitand(crc, 255));
    crc_hi = uint8(bitand(bitshift(crc, -8), 255));
    
    frame = uint8([hex2dec('FD'), header, payload, crc_lo, crc_hi]);
    
    % Increment sequence
    seq = uint8(mod(double(seq) + 1, 256));
end

function crc = mavlink_crc16(data)
% MAVLink X.25 CRC-16 (CRC-16/MCRF4XX)
    crc = uint16(hex2dec('FFFF'));
    for i = 1:length(data)
        tmp = bitxor(uint8(data(i)), uint8(bitand(crc, uint16(255))));
        tmp = bitxor(tmp, uint8(bitand(bitshift(uint16(tmp), 4), uint16(255))));
        t16 = uint16(tmp);
        crc = bitxor( bitshift(crc, -8), ...
              bitxor( bitshift(t16, 8), ...
              bitxor( bitshift(t16, 3), ...
                      bitshift(t16, -4) )));
    end
end

function payload = build_heartbeat_payload()
% HEARTBEAT (msg ID 0): type=6 (GCS), autopilot=8 (invalid),
% base_mode=0, custom_mode=0, system_status=4 (active)
    custom_mode = typecast(uint32(0), 'uint8');     % 4 bytes
    type        = uint8(6);                         % MAV_TYPE_GCS
    autopilot   = uint8(8);                         % MAV_AUTOPILOT_INVALID
    base_mode   = uint8(0);
    sys_status  = uint8(4);                         % MAV_STATE_ACTIVE
    mav_version = uint8(3);                         % MAVLink v2
    payload = uint8([custom_mode, type, autopilot, base_mode, ...
                     sys_status, mav_version]);
end

function payload = build_hil_sensor_payload( ...
        time_usec, xacc, yacc, zacc, xgyro, ygyro, zgyro, ...
        xmag, ymag, zmag, abs_press, diff_press, press_alt, ...
        temperature, fields_updated)
% HIL_SENSOR (msg ID 107) – 64 bytes payload
    payload = uint8(zeros(1, 64));
    payload(1:8)   = typecast(time_usec,      'uint8');
    payload(9:12)  = typecast(xacc,           'uint8');
    payload(13:16) = typecast(yacc,           'uint8');
    payload(17:20) = typecast(zacc,           'uint8');
    payload(21:24) = typecast(xgyro,          'uint8');
    payload(25:28) = typecast(ygyro,          'uint8');
    payload(29:32) = typecast(zgyro,          'uint8');
    payload(33:36) = typecast(xmag,           'uint8');
    payload(37:40) = typecast(ymag,           'uint8');
    payload(41:44) = typecast(zmag,           'uint8');
    payload(45:48) = typecast(abs_press,      'uint8');
    payload(49:52) = typecast(diff_press,     'uint8');
    payload(53:56) = typecast(press_alt,      'uint8');
    payload(57:60) = typecast(temperature,    'uint8');
    payload(61:64) = typecast(fields_updated, 'uint8');
end

function payload = build_hil_gps_payload( ...
        time_usec, lat, lon, alt, eph, epv, vel, ...
        vn, ve, vd, cog, fix_type, satellites)
% HIL_GPS (msg ID 113) – 36 bytes payload
    payload = uint8(zeros(1, 36));
    payload(1:8)   = typecast(time_usec, 'uint8');
    payload(9:12)  = typecast(lat,       'uint8');   % int32 degE7
    payload(13:16) = typecast(lon,       'uint8');
    payload(17:20) = typecast(alt,       'uint8');   % int32 mm
    payload(21:22) = typecast(eph,       'uint8');   % uint16 cm
    payload(23:24) = typecast(epv,       'uint8');
    payload(25:26) = typecast(vel,       'uint8');   % uint16 cm/s
    payload(27:28) = typecast(vn,        'uint8');   % int16 cm/s
    payload(29:30) = typecast(ve,        'uint8');
    payload(31:32) = typecast(vd,        'uint8');
    payload(33:34) = typecast(cog,       'uint8');   % uint16 cdeg
    payload(35)    = fix_type;
    payload(36)    = satellites;
end

function val = bytes_to_single(b4)
% Convert 4 uint8 bytes (little-endian) to a single-precision float.
    val = typecast(uint8(b4), 'single');
end