function q = eul2Quat(eul)
%EUL2QUAT Summary of this function goes here
%   Detailed explanation goes here
% eul = [phi, pitch, yaw]
% ZYX -> Pitch, Yaw, Roll

    eul = [eul(2), eul(3), eul(1)];
    c = cos(eul/2);
    s = sin(eul/2);

    
    % Construct quaternion
    q = [c(:,1).*c(:,2).*c(:,3)+s(:,1).*s(:,2).*s(:,3), ...
          c(:,1).*c(:,2).*s(:,3)-s(:,1).*s(:,2).*c(:,3), ...
          c(:,1).*s(:,2).*c(:,3)+s(:,1).*c(:,2).*s(:,3), ...
          s(:,1).*c(:,2).*c(:,3)-c(:,1).*s(:,2).*s(:,3)];
        
end

