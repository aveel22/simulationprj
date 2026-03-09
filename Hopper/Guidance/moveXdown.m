function ref_new = moveXdown(ref, p, t, tstop, dt)
    % MOVEXDOWN Summary of this function goes here
    %   Detailed explanation goes here
    
    new_params = p;
    new_params.vx = -p.vxl;
    ref_new = moveX(ref, new_params, t, tstop, dt);
end

