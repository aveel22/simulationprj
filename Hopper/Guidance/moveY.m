function ref_new = moveY(ref, p, t, tstop, dt)
    % MOVEX Summary of this function goes here
    %   Detailed explanation goes here
    
    ref_new = ref;
    if t > 0 && t <= tstop
        ref_new(2) = ref(2) + p.vy*dt;
        ref_new(5) = p.vy;
    else
        ref_new(5) = 0;
    end
end

