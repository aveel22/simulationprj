function [hf, hm] = ForceOnComponents(parameters, rc)
%FORCEONCOMPONENTS Summary of this function goes here
%   Detailed explanation goes here
% parameters - geomentry
% rc = [xc, yc, zc] - center of mass

% Single force

eul_f = -deg2rad(parameters(4:6));
eul_L = -deg2rad(parameters(7:9));
hf = Nw(eul_f) * [1;0;0];
L = Nw(eul_L) * [parameters(1); 0; 0];
r = L + [parameters(2); 0; 0];
hm = skew(r - rc') * hf; 