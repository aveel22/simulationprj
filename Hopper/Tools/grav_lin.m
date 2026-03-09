function A = grav_lin(eul, g)
% Возвращает 3x3 матрицу производных ускорения от гравитации
% по углам (phi, theta, psi) в body frame
%
% Параметры:
%   phi   — угол крена (в радианах)
%   theta — угол тангажа (в радианах)
%   psi   — угол рыскания (в радианах)
%   g     — ускорение свободного падения (>0)
%
% Выход:
%   A     — 3x3 матрица: d(ag_body)/d(phi, theta, psi)

phi = eul(1); theta = eul(2); psi = eul(3);

% Предварительные значения
sphi = sin(phi);   cphi = cos(phi);
stheta = sin(theta); ctheta = cos(theta);
spsi = sin(psi);   cpsi = cos(psi);

% Строки по производным
A11 = 0;
A12 = -g * cpsi * ctheta;
A13 =  g * spsi * stheta;

A21 =  g * (sphi * ctheta - spsi * stheta * cphi);
A22 =  g * (-sphi * spsi * ctheta + stheta * cphi);
A23 = -g * sphi * stheta * cpsi;

A31 =  g * (sphi * spsi * stheta + cphi * ctheta);
A32 = -g * (sphi * stheta + spsi * cphi * ctheta);
A33 = -g * stheta * cphi * cpsi;

% Сборка в матрицу
A = [A11, A12, A13;
     A21, A22, A23;
     A31, A32, A33];
end