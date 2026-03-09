function [A,B] = cascade_lsys(n, T)
    % A,B for n cascaded first-order lags with time constant T
    A = (-1/T)*eye(n) + (1/T)*diag(ones(n-1,1), -1);
    B = zeros(n,1); 
    B(1) = 1/T;
end


