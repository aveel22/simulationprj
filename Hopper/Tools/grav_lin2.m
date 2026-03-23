function A = grav_lin2(eul, g)
    phi   = eul(1);
    theta = eul(2);
    psi   = eul(3);

    sphi = sin(phi);   cphi = cos(phi);
    sth  = sin(theta); cth  = cos(theta);
    spsi = sin(psi);   cpsi = cos(psi);

    A11 = 0;
    A12 =  g * sth * cpsi;
    A13 =  g * cth * spsi;

    A21 = -g * (cphi * sth * cpsi + sphi * spsi);
    A22 = -g * (sphi * cth * cpsi);
    A23 =  g * (sphi * sth * spsi + cphi * cpsi);

    A31 =  g * (sphi * sth * cpsi - cphi * spsi);
    A32 = -g * (cphi * cth * cpsi);
    A33 =  g * (cphi * sth * spsi - sphi * cpsi);

    A = [A11, A12, A13;
         A21, A22, A23;
         A31, A32, A33];
end
