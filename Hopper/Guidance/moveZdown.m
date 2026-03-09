function ref_new = moveZdown(ref, p, t, tstop, dt)
    % MOVEXDOWN Summary of this function goes here
    %   Detailed explanation goes here
    
    new_params = p;
    new_params.vz = -p.vz;
    ref_new = moveZ(ref, new_params, t, tstop, dt);
end

