function rot = E2W(ang)
    % E2W Summary of this function goes here
    % Euler's kinematic for NED
    % Yaw(Z) -> Pitch(Y) -> Roll(Z)
    roll = ang(1);
    pitch = ang(2);
    % yaw = ang(3);
    rot = [
        cos(pitch), sin(roll) * sin(pitch),  cos(roll) * sin(pitch);
        0,          cos(roll) * cos(pitch), -sin(roll) * cos(pitch);
        0,          sin(roll),              cos(roll);
        ]  / cos(pitch);
end

