function new_params = disturb_parameters(params, w)
    %DISTURB_PARAMETERS Summary of this function goes here
    %   Detailed explanation goes here
    % Copy data
    new_params = params; 
    % Disturb data
    new_params.Inertia = params.Inertia * (1 + w * (2*rand - 1));
    new_params.mass = params.mass * (1 + w * (2*rand - 1));
    new_params.rcgs = params.rcgs * (1 + w * (2*rand - 1));
    new_params.Jrot = params.Jrot * (1 + w * (2*rand - 1));
    new_params.Cv = params.Cv * (1 + w * (2*rand - 1));
    new_params.Cw = params.Cw * (1 + w * (2*rand - 1));
    new_params.Jdot_mlt = params.Jdot_mlt * (1 + w * (2*rand - 1));
    new_params.m = params.m * (1 + w * (2*rand - 1));
    new_params.beta = params.beta * (1 + w * (2*rand - 1));
    new_params.gamma = params.gamma * (1 + w * (2*rand - 1));
    new_params.Jdot = params.Jdot * (1 + w * (2*rand - 1));
    
    new_params.Jo = new_params.Inertia(:,:,1);
    new_params.geom = distirb_geom(params.geom);
    [Hf, Hm] = ComponentForForces(geometry, params.rcgs(1,:));
    new_params.Hf = Hf;
    new_params.Hm = Hm;
end

