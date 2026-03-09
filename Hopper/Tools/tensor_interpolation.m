function T_interp = tensor_interpolation(t_input, mass, Tensors)
    % Load the tensor data

    [N, M, ~] = size(Tensors); % N x M matrices over Count time steps

    % Preallocate the interpolated tensor
    T_interp = zeros(N, M);

    % Interpolate each tensor component
    for i = 1:N
        for j = 1:M
            T_interp(i,j) = interp1(mass, squeeze(Tensors(i, j, :)), t_input, 'spline');
        end
    end
end


