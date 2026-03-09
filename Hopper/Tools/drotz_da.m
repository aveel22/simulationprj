function dR = drotz_da(a)
    % Derivative of rotz w.r.t. angle a
    c = cos(a);
    s = sin(a);
    dR = [-s,  c,  0;
          -c, -s,  0;
           0,  0,  0];
end