function vector = Tensor2Vector(tensor)
%TENSOR2VECTOR Summary of this function goes here
%   Detailed explanation goes here
vector = [tensor(1,1), tensor(2,2), tensor(3,3), ...
    tensor(1,2), tensor(1,3), tensor(2,3)];

