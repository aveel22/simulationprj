function TensorDiag = DiagTensorInertia(J)
%% make diagonal tensor of inertia
% J - tensor of inertia

if size(J,1) == 1 | size(J,2) == 1
    d = diag(J);
elseif size(J,1) == 3 & size(J,2) == 3
    d = diag(diag(J));
else 
    d = J;
end
TensorDiag = d;