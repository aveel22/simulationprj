%% Start Script
rxd = 0;     ryd = 0;       rzd = 0;
vxd = 0;     vyd = 0;       vzd = 0;
phid = 0;   thetad = pi/2;  psid = 0;
wxd = 0;    wyd = 0;        wzd = 0;

roll_max = 180;
roll_wx_max = 5;

dtau = 0.1;
t_pause = 5;
ref1 = [rxd; ryd; rzd; vxd; vyd; vzd; phid; thetad; psid; wxd; wyd; wzd];
dim = size(ref1, 1);

%% Define the desired trajectory as a function of time
%% 1. Vertical ascending. Tf = 25
Hmax = 50;
vxd1 = 0;
vyd1 = 2;
vzd1 = 0;
Tf1 = Hmax / vyd1;
Tspan = 0:dtau:Tf1;
desired_trajectory = @(t) [rxd; vyd1*t; rzd; 
                            vxd1; vyd1; vzd1; 
                            phid; thetad; psid; 
                            wxd; wyd; wzd];

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory(Tspan(i));
end

% The 1st part of trajectory
sim_1 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end)+t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(5) = 0;
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_1 = [sim_1; temp(2:end,:)];

%% 2. Shifting in positive X direction
Tf2 = 10;
vxd2 = 0.5;
vyd2 = 0;
vzd2 = 0;
Tspan = dtau:dtau:Tf2;
desired_trajectory2 = @(t) [vxd2*t; vyd1*Tf1; rzd; 
                            vxd2; vyd2; vzd2; 
                            phid; thetad; psid; 
                            wxd; wyd; wzd];

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory2(Tspan(i));
end

Tspan = Tspan + sim_1(end,1);
% The 2nd part of trajectory
sim_2 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(4) = 0;                   % make vxd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_2 = [sim_2; temp(2:end,:)];

%% 3. Shifting in negative Z direction
Tf3 = 10;
vxd3 = 0;
vyd3 = 0;
vzd3 = -0.5;
Tspan = dtau:dtau:Tf3;
desired_trajectory3 = @(t) [vxd2*Tf2; vyd1*Tf1; vzd3*t; 
                            vxd3; vyd3; vzd3; 
                            phid; thetad; psid; 
                            wxd; wyd; wzd];
dim = size(ref1,1);

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory3(Tspan(i));
end

Tspan = Tspan + sim_2(end,1);
% The 3rd part of trajectory
sim_3 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(6) = 0;                   % make vzd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_3 = [sim_3; temp(2:end,:)];

%% 4. Shifting back in negative X direction
Tf4 = 10;
vxd4 = -0.5;
vyd4 = 0;
vzd4 = 0;
Tspan = dtau:dtau:Tf4;
desired_trajectory4 = @(t) [vxd2*Tf2+vxd4*t; vyd1*Tf1; vzd3*Tf3; 
                            vxd4; vyd4; vzd4; 
                            phid; thetad; psid; 
                            wxd; wyd; wzd];

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory4(Tspan(i));
end

Tspan = Tspan + sim_3(end,1);
% The 4th part of trajectory
sim_4 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(4) = 0;                   % make vxd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_4 = [sim_4; temp(2:end,:)];

%% 5. Shifting back in positive Z direction
Tf5 = 10;
vxd5 = 0;
vyd5 = 0;
vzd5 = 0.5;
Tspan = dtau:dtau:Tf5;
desired_trajectory5 = @(t) [vxd2*Tf2+vxd4*Tf4; vyd1*Tf1; vzd3*Tf3+vzd5*t; 
                            vxd5; vyd5; vzd5; 
                            phid; thetad; psid; 
                            wxd; wyd; wzd];

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory5(Tspan(i));
end

Tspan = Tspan + sim_4(end,1);
% The 5th part of trajectory
sim_5 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(6) = 0;                   % make vzd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_5 = [sim_5; temp(2:end,:)];

%% 6. Roll around 
vxd6 = 0;
vyd6 = 0;
vzd6 = 0;
phi_max = deg2rad(roll_max);
wxd6 = deg2rad(roll_wx_max);
Tf6 = abs(phi_max / wxd6);
Tspan = dtau:dtau:Tf6;
                         
desired_trajectory6 = @(t) [vxd2*Tf2+vxd4*Tf4; vyd1*Tf1; vzd3*Tf3+vzd5*Tf5; 
                            vxd6; vyd6; vzd6; 
                            wxd6*t; thetad; psid; 
                            wxd6; wyd; wzd];

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory6(Tspan(i));
end

Tspan = Tspan + sim_5(end,1);
% The 6th part of trajectory
sim_6 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(10) = 0;                   % make wxd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_6 = [sim_6; temp(2:end,:)];

%% 7. Rolling back
vxd7 = 0;
vyd7 = 0;
vzd7 = 0;
phi_max = deg2rad(roll_max);
wxd7 = deg2rad(-roll_wx_max);
Tf7 = abs(phi_max / wxd7);

Tspan = dtau:dtau:Tf7;
desired_trajectory7 = @(t) [vxd2*Tf2+vxd4*Tf4; vyd1*Tf1; vzd3*Tf3+vzd5*Tf5; 
                            vxd7; vyd7; vzd7; 
                            wxd6 * Tf6 + wxd7 * t; thetad; psid; 
                            wxd7; wyd; wzd];

trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory7(Tspan(i));
end

Tspan = Tspan + sim_6(end,1);
% The 7th part of trajectory
sim_7 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause;
trj_pause = trajectory_1(end,:);
trj_pause(10) = 0;                   % make wxd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_7 = [sim_7; temp(2:end,:)];

%% 8. Descending (in negative Y direction)
vxd8 = 0;
vyd8 = -2;
vzd8 = 0;
Tf8 = 20;
Tspan = dtau:dtau:Tf8;

desired_trajectory8 = @(t) [vxd2*Tf2+vxd4*Tf4; vyd1*Tf1+vyd8*t; vzd3*Tf3+vzd5*Tf5; 
                            vxd8; vyd8; vzd8; 
                            wxd6 * Tf6 + wxd7 * Tf7; thetad; psid; 
                            wxd; wyd; wzd];

% dim = size(ref1, 1);
trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory8(Tspan(i));
end

Tspan = Tspan + sim_7(end,1);
% The 8th part of trajectory
sim_8 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause / 5;
trj_pause = trajectory_1(end,:);
trj_pause(5) = 0;                   % make vzd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_8 = [sim_8; temp(2:end,:)];

%% 9. Landing (in negative Y direction)
vxd9 = 0;
vyd9 = -0.5;
vzd9 = 0;
Tf9 = abs(sim_8(end, 3) / vyd9);
Tspan = dtau:dtau:Tf9;

desired_trajectory9 = @(t) [vxd2*Tf2+vxd4*Tf4; vyd1*Tf1+vyd8*Tf8 + vyd9*t; vzd3*Tf3+vzd5*Tf5; 
                            vxd9; vyd9; vzd9; 
                            wxd6 * Tf6 + wxd7 * Tf7; thetad; psid; 
                            wxd; wyd; wzd];
% dim = size(ref1, 1);
trajectory_1 = zeros(length(Tspan), dim); % Preallocate for efficiency
for i = 1:length(Tspan)
    trajectory_1(i, :) = desired_trajectory9(Tspan(i));
end

Tspan = Tspan + sim_8(end,1);
% The 9th part of trajectory
sim_9 = [Tspan', trajectory_1];
T1 = Tspan(end):dtau:Tspan(end) + t_pause / 5;
trj_pause = trajectory_1(end,:);
trj_pause(5) = 0;                   % make vzd = 0
temp = [T1', repmat(trj_pause, length(T1), 1)];
sim_9 = [sim_9; temp(2:end,:)];

%% Reference Trajectory
% delete trajectory_1
sim = [sim_1; sim_2; sim_3];
ref_traj = [sim_1; sim_2; sim_3; sim_4; sim_5; sim_6; sim_7; sim_8; sim_9];
%% Extra
% X = sim(:, [1, 2]);
% Y = sim(:, [1, 3]);
% Z = sim(:, [1, 4]);
% Roll = sim(:, [1, 8]);
% Pitch = sim(:, [1, 9]);
% Yaw = sim(:, [1, 10]);
%% Delete obsolete data
clear ref1 rxd ryd rzd vxd vyd vzd phid thetad psid wxd wyd wzd;
clear roll_max roll_wx_max dtau = 0.1 t_pause dim;
clear Hmax vxd1 vyd1 vzd1 Tf1 desired_trajectory;
clear Tf2 vxd2 vyd2 vzd2 desired_trajectory2;
clear Tf3 vxd3 vyd3 vzd3 desired_trajectory3;
clear Tf4 vxd4 vyd4 vzd4 desired_trajectory4;
clear Tf5 vxd5 vyd5 vzd5 desired_trajectory5;
clear vxd6 vyd6 vzd6 phi_max wxd6 Tf6 desired_trajectory6;
clear vxd7 vyd7 vzd7 phi_max wxd7 Tf7 desired_trajectory7;
clear vxd8 vyd8 vzd8 Tf8 desired_trajectory8;
clear vxd9 vyd9 vzd9 Tf9 Tspan desired_trajectory9;
clear Tspan T1 trj_pause trajectory_1 temp
clear sim_1 sim_2 sim_3 sim_4 sim_5 sim_6 sim_7 sim_8 sim_9 i

