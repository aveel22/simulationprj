function y = rotx(a)
    % ROTX Build rotation matrix around X axis
    % x - vector   
    y = [1,     0,      0; 
         0, cos(a), sin(a); 
         0, -sin(a), cos(a)];
