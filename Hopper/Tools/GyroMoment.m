function L = GyroMoment(omega, xi, params)
    %% return disturbance force and moments caused by gyromoments
    % omega = [wx; wy; wz] - body rotation vector
    % xi = [xi_x, xi_y, xi_z] - rotor angular velocity
    % Jrot - moment inertia of jet's rotor
    L = zeros(3,1);
    Jrot = params.Jrot;
    %% Kinetic moment
    %  Krot = Jrot * xi;
    L(1:3, 1) = cross(omega, Jrot * xi);