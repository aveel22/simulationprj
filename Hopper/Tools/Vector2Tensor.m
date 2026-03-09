function tensor = Vector2Tensor(vector)
%VECTOR2TENSOR Summary of this function goes here
%   Detailed explanation goes here
tensor = [vector(1), vector(4), vector(5);
          vector(4), vector(2), vector(6);
          vector(5), vector(6), vector(3)];

