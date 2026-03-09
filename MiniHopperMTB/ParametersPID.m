% Attitude PID
pid_params_att = struct();
u_mlt = 1;
dt = 0.1;
pid_params_att.p = u_mlt .* [1.; 1.; 0.5] * 2.;
pid_params_att.i = u_mlt .* [1.; 1; 1.] * 1e-2;
pid_params_att.d = u_mlt .* [1.; 1.; 0.5] * 2.;
pid_params_att.ff = u_mlt .* [1.; 1.; 0.5] * 1;
pid_params_att.dt = dt;

% Rates PID
pid_params_rate = struct();
pid_params_rate.p = u_mlt .* [1.; 1.; 0.5] * 2.;
pid_params_rate.i = u_mlt .* [1.; 1; 1.] * 1e-2;
pid_params_rate.d = u_mlt .* [1.; 1.; 0.5] * 2.;
pid_params_rate.ff = u_mlt .* [1.; 1.; 0.5] * 1;
pid_params_rate.dt = dt;

% Position PID
pid_params_pos = struct();
pid_params_pos.p = u_mlt .* [1.; 1.; 1] * 2.;
pid_params_pos.i = u_mlt .* [1.; 1; 10.] * 1e-2;
pid_params_pos.d = u_mlt .* [1.; 1.; 0.5] * 2.;
pid_params_pos.i_a = 2;
pid_params_pos.i_acc = 5;
pid_params_pos.dt = dt;

% Altitude PID
pid_params_alt = struct();
pid_params_alt.p = u_mlt .* [1.; 1.; 1] * 2.;
pid_params_alt.i = u_mlt .* [1.; 1; 10.] * 1e-2;
pid_params_alt.d = u_mlt .* [1.; 1.; 0.5] * 2.;
pid_params_alt.i_a = 2;
pid_params_alt.i_acc = 5;
pid_params_alt.dt = dt;