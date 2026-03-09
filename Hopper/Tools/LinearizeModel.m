function [A, B, C, D] = LinearizeModel(x,u,p)
%LINEARIZEMODEL Summary of this function goes here
%   Detailed explanation goes here
    [A, B] = build_lsys(x, u, p);
    C = eye(size(A,1));
    D = zeros(size(A, 1), size(B, 2));
end

