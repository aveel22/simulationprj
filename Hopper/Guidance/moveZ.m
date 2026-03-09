function ref_new = moveZ(ref, p, t, tstop, dt)
    % MOVEX Summary of this function goes here
    %   Detailed explanation goes here
    
    ref_new = ref;
    if t > 0 && t < tstop
        ref_new(3) = ref(3) + p.vz*dt;
        ref_new(6) = p.vz;
    else
        ref_new(6) = 0;
    end
end

