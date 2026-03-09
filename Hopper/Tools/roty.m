function y = roty(a)
    % ROTY Buiold rotation matrix around Y axis
    % x - vector   
    y = [cos(a),    0,  -sin(a); 
         0,         1,       0; 
         sin(a),    0,   cos(a)];