function res = skew(arg)
    x = arg(1);
    y = arg(2);
    z = arg(3);
    res = [0, -z,  y;
           z,  0, -x;
           -y, x,  0];
end

