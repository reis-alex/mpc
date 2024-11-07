function [P0_d] = P0_dot(R0,t0,robot)
%P0_DOT Summary of this function goes here
%   Detailed explanation goes here
    P0_d = [SkewSym(t0(1:3))*R0  zeros(3);
            zeros(3)             zeros(3)];
end