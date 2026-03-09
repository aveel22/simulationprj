function [Hf, Hm] = ComponentForForces(geom, rc, skip_column, Nrot);
%COMPONENTFORFORCES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    Nrot = eye(3);
end
if nargin < 3
    skip_column = [];
end

Hf = [];
Hm = [];
k = 1;
for i=1:9
    if ismember(i,skip_column)
        continue;
    end
    if i == 1
        [hf, hm] = JetForceOnComponents(geom(i,:), rc);
    else
        [hf, hm] = ForceOnComponents(geom(i,:), rc);
    end
    hf = Nrot * hf;
    hm = Nrot * hm;
    Hf = horzcat(Hf, hf);
    Hm = horzcat(Hm, hm);
    k = k + 1;
end