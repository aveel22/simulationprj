function [Fxyz_total, Mxyz_total] = ground_contact(Fxyz_motors, Mxyz_motors, Xe, Ve_body, omega_b, DCM_be)
% GROUND_CONTACT  Simple spring-damper ground reaction model
%
% Add as a MATLAB Function block between the throttle-to-force mapping
% and the 6-DOF block.  It adds ground reaction forces/torques when
% the vehicle is at or below ground level.
%
% INPUTS:
%   Fxyz_motors  [3x1]  Motor forces in body frame (FRD) [N]
%   Mxyz_motors  [3x1]  Motor torques in body frame [N·m]
%   Xe           [3x1]  Position in NED [m] — from 6-DOF
%   Ve_body      [3x1]  Body-frame velocity [m/s] — from 6-DOF Vb output
%   omega_b      [3x1]  Body angular rates [rad/s]
%   DCM_be       [3x3]  Body-to-Earth DCM — from 6-DOF
%
% OUTPUTS:
%   Fxyz_total   [3x1]  Forces with ground reaction added (body frame)
%   Mxyz_total   [3x1]  Torques with ground damping added (body frame)

% Ground plane parameters
K_spring  = 5000;    % Normal spring stiffness [N/m]
C_damp    = 500;     % Normal damping [N·s/m]
K_lateral = 2000;    % Lateral friction spring [N/m]
C_lateral = 200;     % Lateral friction damping [N·s/m]
C_rot     = 5.0;     % Rotational damping [N·m·s/rad]

% Start with motor forces/torques
Fxyz_total = Fxyz_motors;
Mxyz_total = Mxyz_motors;

% NED: Z positive = down.  Ground is at Z = 0.
% Penetration depth: positive when below ground
penetration = Xe(3);  % > 0 means below ground

if penetration > 0
    % Earth-frame velocity from body velocity
    Ve_earth = DCM_be' * Ve_body;  % DCM_be is body-to-earth, transpose = earth-to-body inverse
    % Actually DCM_be from Simulink 6-DOF is body-to-earth rotation
    % Ve_earth = DCM_be * Ve_body would be wrong direction
    % The 6-DOF block outputs DCM_be such that: V_earth = DCM_be' * V_body
    % Let me use it correctly:
    Ve_earth = DCM_be' * Ve_body;
    
    % ---- Normal force (Z-axis, push UP = negative in NED) ----
    Fz_ground_earth = -K_spring * penetration - C_damp * max(Ve_earth(3), 0);
    
    % ---- Lateral friction (X,Y axes in earth frame) ----
    Fx_friction_earth = -K_lateral * 0 - C_lateral * Ve_earth(1);
    Fy_friction_earth = -K_lateral * 0 - C_lateral * Ve_earth(2);
    
    % Combined ground force in earth frame
    F_ground_earth = [Fx_friction_earth; Fy_friction_earth; Fz_ground_earth];
    
    % Transform to body frame
    F_ground_body = DCM_be * F_ground_earth;
    
    % ---- Rotational damping (prevents spinning on ground) ----
    M_ground_body = -C_rot * omega_b;
    
    % Add to motor forces/torques
    Fxyz_total = Fxyz_motors + F_ground_body;
    Mxyz_total = Mxyz_motors + M_ground_body;
end

end