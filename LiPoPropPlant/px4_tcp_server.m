function dataout = px4_tcp_server(datain)
    % Persistent handles to manage stream state machines
    persistent server_socket client_socket in_stream out_stream init_done command_sent
    
    % Initialize a flat 512-byte buffer stream matching your MAVLink block size
    dataout = zeros(1, 512, 'double');
    
    % Step 1: First-simulation-step variable initialization
    if isempty(init_done)
        server_socket = []; client_socket = []; in_stream = []; out_stream = [];
        init_done = 1; command_sent = 0;
    end
    
    % Step 2: Spin up the Java TCP Server socket if it is closed
    if isempty(server_socket)
        try
            server_socket = java.net.ServerSocket(4560);
            server_socket.setSoTimeout(10); % Quick 10ms timeout loop window
            
            client_socket = server_socket.accept();
            in_stream = client_socket.getInputStream();
            out_stream = client_socket.getOutputStream();
            disp('🚀 PX4 Autopilot connected successfully to Simulink!');
            command_sent = 0; % Reset automation trigger flag
        catch
            if ~isempty(server_socket)
                try server_socket.close(); catch; end
                server_socket = [];
            end
            return;
        end
    end
    
    % Make sure your Step 3 code looks clean and doesn't inject text strings:
    if ~isempty(out_stream) && ~isempty(datain)
        try
            java_send_buffer = int8(datain);
            out_stream.write(java_send_buffer);
            out_stream.flush();
        catch
            server_socket = []; client_socket = []; in_stream = []; out_stream = [];
            return;
        end
    end

    
    % Step 4: Write standard simulation plant feedback data TO PX4
    if ~isempty(out_stream) && ~isempty(datain) && (command_sent == 1)
        try
            java_send_buffer = int8(datain);
            out_stream.write(java_send_buffer);
            out_stream.flush();
        catch
            server_socket = []; client_socket = []; in_stream = []; out_stream = []; command_sent = 0;
            return;
        end
    end
    
    % Step 5: Read incoming motor control signals FROM PX4
    if ~isempty(in_stream)
        try
            avail = in_stream.available();
            if avail > 0
                bytes_to_read = min(avail, 512);
                buffer = jaggedArray(java.lang.Byte.TYPE, bytes_to_read);
                in_stream.read(buffer);
                
                for i = 1:bytes_to_read
                    b = buffer(i);
                    if b < 0
                        dataout(i) = double(b + 256);
                    else
                        dataout(i) = double(b);
                    end
                end
            end
        catch
            server_socket = []; client_socket = []; in_stream = []; out_stream = []; command_sent = 0;
        end
    end
end

% Quick helper function to build native Java array allocations inside MATLAB
function arr = jaggedArray(type, len)
    arr = java.lang.reflect.Array.newInstance(type, len);
end
