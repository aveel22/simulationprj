% PX4_SIM_CLEANUP  Force-release port 4560 after a crashed simulation
%
% When Simulink crashes or you hit Stop without a clean shutdown, the
% Java ServerSocket may still hold port 4560.  Run this script to kill
% it, then restart your simulation.
%
% Usage:  px4_sim_cleanup

fprintf('Clearing persistent variables in px4_sim_bridge...\n');
clear px4_sim_bridge

% Force Java garbage collection to release the socket
java.lang.System.gc();
pause(0.5);

fprintf('Port 4560 should be free now.\n');
fprintf('If it is still busy, close and reopen MATLAB.\n');