function Fd = GyroDisturb(omega, xi, params)
    %% return disturbance force and moments caused by gyromoments
    % omega = [wx; wy; wz] - body rotation vector
    % xi = [xi_x, xi_y, xi_z] - rotor angular velocity
    Fd = zeros(6,1);
    %% Kinetic moment
    %  Krot = Jrot * xi;
    Fd(4:6, 1) = GyroMoment(omega, xi, params);