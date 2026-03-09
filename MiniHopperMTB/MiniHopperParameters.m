% System parameters
nx = 12; nu = 6;
% Force levers
Lz = 0.2;       % Roll and Pitch lever
Lr = 0.03;       % Yaw (heading) lever

hopper_ode = struct();
hopper_ode.rc = [0; 0; 0.2];
hopper_ode.J = [3, 0, 0; 0, 3, 0; 0, 0, 0.45] * 0.01;
hopper_ode.m = 1.5;
hopper_ode.Cv = eye(3) * 0.1;
hopper_ode.Cw = eye(3) * 0.01;
hopper_ode.Jdot = zeros(3);
hopper_ode.g = 9.8065;
hopper_ode.n = 2;
hopper_ode.T = 0.5;
hopper_ode.gc = 0;
hopper_ode.Lz = Lz;
hopper_ode.Lr = Lr;
% Transform from servo angle to moment (Force = 1 [N])
hopper_ode.D2M = [[Lr, Lr, Lr, Lr];[-Lz, 0, Lz, 0];[0, -Lz, 0, Lz]];
% Transform back from moment to angle position
hopper_ode.M2D = [
                    [0.25 / Lr, -0.5 / Lz,          0]; 
                    [0.25 / Lr,     0,      -0.5 / Lz];
                    [0.25 / Lr, 0.5 / Lz,           0]; 
                    [0.25 / Lr,     0,       0.5 / Lz];
                  ];

hopper_ode.inv.roll = true;
hopper_ode.inv.yaw = true;
hopper_ode.inv.pitch = false;
