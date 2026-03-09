function ref_new = moveX(ref, p, t, tstop, dt)
    % MOVEX Summary of this function goes here
    %   Detailed explanation goes here
    
    ref_new = ref;
    if t > 0 && t <= tstop
        ref_new(1) = ref(1) + p.vx*dt;
        ref_new(4) = p.vx;
    else
        ref_new(4) = 0;
    end
end

