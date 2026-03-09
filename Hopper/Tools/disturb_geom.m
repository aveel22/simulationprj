function new_geom = disturb_geom(geom, w)
    %DISTURB_GEOM Summary of this function goes here
    %   Detailed explanation goes here
    if nargin<2
        w = 0.01;
    end
    new_geom = geom;
    [nrow, ncols] = size(geom);
    for i=1:ncols
        new_geom(:,i) = geom(:,i).*(1+w*(2*rand(nrow,1) - 1));
    end
end

