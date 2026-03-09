function [x_all, y_all, y1_all] = merge_with_interpolation(x, y, x1, y1)
    % Merge two sets of (x, y) points:
    % - (x, y)  : base grid and values
    % - (x1, y1): extra grid and values
    %
    % Returns:
    % x_all  : merged sorted x values
    % y_all  : original y extended with interpolation at x1
    % y1_all : original y1 extended with zeros at x
    
    % Step 1: merge and sort unique x values
    x_all = unique([x(:); x1(:)]);  % ensure column vectors, merge & sort
    
    % Step 2: interpolate y onto the merged grid
    y_all = interp1(x, y, x_all, 'linear');
    
    % Step 3: extend y1 with zeros at original x locations
    y1_all = zeros(size(x_all));
    
    % Find positions of x1 inside x_all and assign y1 values
    [~, loc] = ismember(x1, x_all);
    y1_all(loc) = y1;
end
