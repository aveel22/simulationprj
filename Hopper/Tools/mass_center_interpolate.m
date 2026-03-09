function rcg_interp = mass_center_interpolate(m_input, mass, rcgs)
%% Center of mass interpolation
% m_input   - current mass
% mass      - array of masses
% rcgs      - array of center of mass according to `mass` array

    % Preallocate the interpolated center of mass array
    rcg_interp = zeros(1, 3); % One row with x, y, z

    % Interpolate each coordinate (x, y, z) separately
    for j = 1:3
        rcg_interp(j) = interp1(mass, rcgs(:, j), m_input, 'spline');
    end
end


