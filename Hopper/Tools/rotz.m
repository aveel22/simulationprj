function y = rotz(a)
    % ROTZ Build rotation matrix around Z
    % x - vector   
    y = [cos(a),  sin(a),   0; 
         -sin(a), cos(a),   0; 
         0,           0,    1];
