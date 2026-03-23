function R = Nw2(eul)
    phi   = eul(1); % roll
    theta = eul(2); % pitch
    psi   = eul(3); % yaw

    cphi = cos(phi);   sphi = sin(phi);
    cth  = cos(theta); sth  = sin(theta);
    cpsi = cos(psi);   spsi = sin(psi);

    R = [ ...
        cth*cpsi,                    cth*spsi,                   -sth;
        sphi*sth*cpsi - cphi*spsi,   sphi*sth*spsi + cphi*cpsi,  sphi*cth;
        cphi*sth*cpsi + sphi*spsi,   cphi*sth*spsi - sphi*cpsi,  cphi*cth];
end

