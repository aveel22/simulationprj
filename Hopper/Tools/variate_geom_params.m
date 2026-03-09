function newparams = variate_geom_params(param, w)
    %VARIATE_GEOM_PARAMS Summary of this function goes here
    %   Detailed explanation goes here
    
    newparams = param;
    dim = size(param.geom);
    for i = 1:dim(2)
        f = w(i) * (2*rand(dim(1),1) - 1);
        newparams.geom(:,i) = newparams.geom(:,i) + f;
    end
end
