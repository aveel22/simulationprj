function ref_new = moveYdown(ref, p, t, tstop, dt)
    % MOVEXDOWN Summary of this function goes here
    %   Detailed explanation goes here
    
    new_params = p;
    new_params.vy = -p.vy;
    ref_new = moveY(ref, new_params, t, tstop, dt);
end

