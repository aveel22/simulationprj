function J_rotated = RotateTensor(J, eul)
%% Rotation tensor of inertia
% J     - tensor of inertia
% eul   - Euler's angles

%% Unpack Euler's angles
phi= eul(1);
theta = eul(2);
psi = eul(3);
%% Build Rotation Matrix
R = rotz(theta) * roty(psi) * rotx(phi);
%% Tensor Rotation
J_rotated = R * J * R';