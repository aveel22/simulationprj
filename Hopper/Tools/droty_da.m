function dR = droty_da(a)
    % Derivative of roty w.r.t. angle a
    c = cos(a);
    s = sin(a);
    dR = [-s,  0, -c;
           0,  0,  0;
           c,  0, -s];
end

