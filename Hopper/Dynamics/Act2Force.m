function F = Act2Force(u, p)
% Transform actuators forces to 6DoF Forces
% u - Control Forces: F1, F2, F3, F4, F5, F6, F7, F8, F9
% p - parameters
%% Unpack data
F = zeros(6,1);
h = p(1);
%%
Hf = [1, 1, 1, 1, 0, 0, 0, 0, 1;
     -1, 0, 1, 0, 0, 0, 0, 0, 0;
      0,-1, 0, 1, 0, 0, 0, 0, 0];
% Lever matrix
%     F1 F2 F3 F4 F5 F6 F7 F8 F9
Hm = [h, h, h, h, 0, 0, 0,  0, 0;
      0, 0, 0, 0, h, 0, -h, 0, 0;
      0, 0, 0, 0, 0, h, 0, -h, 0];
F(1:3) = Hf * u;
F(4:6) = Hm * u;