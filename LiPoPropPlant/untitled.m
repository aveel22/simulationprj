%% thrust_linearization.m
% Physics-based throttle-to-thrust linearization for Simulink
% Motor: Brother Hobby T5 3115-1200KV + Gemfan Vortex 7x7 + 6S
%
% Implements two approaches:
%   1. Physics-based power law:  delta = T_norm^n   (n=0.5613 from params)
%   2. PX4-style with voltage compensation (best for full simulation)
%
% Usage in Simulink:
%   - Use as MATLAB Function block
%   - Input:  T_norm  (normalized thrust 0→1)
%             U_bat   (current battery voltage, V)
%   - Output: delta   (ESC duty cycle 0→1)

function delta = thrust_to_delta(T_norm, U_bat)
    % ── Identified parameters ──────────────────────────────────
    K_e   = 0.007958;   % V·s/rad
    C_t   = 4.70e-6;    % N·s²/rad²
    C_d   = 5.19e-8;    % N·m·s²/rad²
    K_t   = K_e;        % BLDC law
    R_mot = 0.127;      % Ω
    R_esc = 0.047;      % Ω
    R_tot = R_mot + R_esc;
    U_nom = 22.2;       % 6S nominal voltage (V)

    % ── Max thrust at current battery voltage ──────────────────
    % At delta=1, steady-state:  U_bat = K_e·ω + (C_d·ω²/K_t)·R_tot
    % Solve quadratic for ω_max:  (C_d·R_tot/K_t)·ω² + K_e·ω - U_bat = 0
    a = C_d * R_tot / K_t;
    b = K_e;
    c = -U_bat;  % uses CURRENT battery voltage, not nominal
    omega_max = (-b + sqrt(max(b^2 - 4*a*c, 0))) / (2*a);
    T_max = C_t * omega_max^2;

    % ── Desired absolute thrust ────────────────────────────────
    T_desired = T_norm * T_max;

    % ── Invert: find omega that produces T_desired ─────────────
    % T = C_t·ω²  →  ω_desired = √(T/C_t)
    omega_desired = sqrt(max(T_desired, 0) / C_t);

    % ── Invert motor model: find delta from omega ──────────────
    % δ·U_bat = K_e·ω + I·R_tot  where I = C_d·ω²/K_t
    V_required = K_e * omega_desired + (C_d * omega_desired^2 / K_t) * R_tot;
    delta = V_required / max(U_bat, 1.0);  % clamp U_bat>0

    % ── Saturate to valid range ────────────────────────────────
    delta = max(0.0, min(1.0, delta));
end

%% ── Simplified version (no U_bat input, for initial testing) ──
% Uses pre-fitted power law: delta = T_norm^0.5613
% Good approximation at nominal voltage, no battery feedback needed

function delta = thrust_to_delta_simple(T_norm)
    n     = 0.5613;   % fitted exponent from identified parameters
    delta = max(0.0, min(1.0, T_norm^n));
end

%% ── Lookup table version (most Simulink-friendly) ─────────────
% Pre-compute at nominal voltage for use as 1D lookup table in Simulink
% Copy T_norm_lut and delta_lut vectors into a 1D Lookup Table block

function [T_norm_lut, delta_lut] = compute_lut(U_bat_nominal)
    if nargin < 1, U_bat_nominal = 22.2; end

    K_e   = 0.007958; C_t = 4.70e-6; C_d = 5.19e-8;
    K_t   = K_e; R_mot = 0.127; R_esc = 0.047;
    R_tot = R_mot + R_esc; U_bat = U_bat_nominal;

    % Sweep delta 0→1, compute T_norm at each point
    deltas = linspace(0, 1, 200);
    T_norm_lut = zeros(size(deltas));
    a = C_d * R_tot / K_t;
    b = K_e;
    for i = 1:length(deltas)
        d = deltas(i);
        if d < 1e-4; T_norm_lut(i) = 0; continue; end
        c_coef = -d * U_bat;
        omega  = (-b + sqrt(max(b^2 - 4*a*c_coef, 0))) / (2*a);
        T_norm_lut(i) = 4.70e-6 * omega^2;
    end
    T_max      = T_norm_lut(end);
    T_norm_lut = T_norm_lut / T_max;
    delta_lut  = deltas;

    fprintf('Lookup table computed at U_bat = %.1f V\n', U_bat_nominal);
    fprintf('T_max (single motor) = %.3f N\n', T_max);
    fprintf('Copy T_norm_lut (breakpoints) and delta_lut (table data)\n');
    fprintf('into a 1D Lookup Table block in Simulink.\n');
end