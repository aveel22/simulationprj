function dR = drotx_da(a)
    % Derivative of rotx w.r.t. angle a
    c = cos(a);
    s = sin(a);
    dR = [0,  0,  0;
          0, -s,  c;
          0, -c, -s];
end

