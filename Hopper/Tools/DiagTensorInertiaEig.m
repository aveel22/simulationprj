function TensorDiag = DiagTensorInertiaEig(J)
    %% Diagonalize Tensor of Inertia
    % J - inertia tensor (should be 3x3)

    % Check if J is 3x3
    if size(J,1) == 3 && size(J,2) == 3
        % Perform eigenvalue decomposition to diagonalize
        [V, D] = eig(J); 
        TensorDiag = D;  % D contains the diagonalized values
    else
        error('Input must be a 3x3 matrix');
    end
end


