function [hf, hm] = JetForceOnComponents(parameters, rc)
%JETFORCEONCOMPONENTS Summary of this function goes here
%   Detailed explanation goes here

eul_f = -deg2rad(parameters(4:6));
eul_L = -deg2rad(parameters(7:9));
hf = Nw(eul_f) * [1;0;0];
L = Nw(eul_L) * [0; parameters(2); parameters(3)];
r = L + [parameters(1); 0; 0];
hm = skew(r - rc') * hf;