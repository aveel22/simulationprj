function T = Tw2(eul)
    phi   = eul(1); % roll
    theta = eul(2); % pitch

    ctheta = cos(theta);
    if abs(ctheta) < 1e-3
        ctheta = 1e-3 * (2*(ctheta >= 0)-1);
    end

    T = [ ...
        1,  sin(phi)*tan(theta),  cos(phi)*tan(theta);
        0,  cos(phi),            -sin(phi);
        0,  sin(phi)/ctheta,      cos(phi)/ctheta ];
end

