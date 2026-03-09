function M = QM(q)
%QM Summary of this function goes here
%   Detailed explanation goes here

q_0 = q(1);
q_1 = q(2);
q_2 = q(3);
q_3 = q(4);

M = [  -q_1, -q_2, -q_3;
        q_0, -q_3, q_2;
        q_3, q_0, -q_1;
        -q_2, q_1, q_0];
