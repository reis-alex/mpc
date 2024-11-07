function P0_d0 = P0_dr0(R0)
%P0_DR0 Summary of this function goes here
%   Detailed explanation goes here
P0_d0 = zeros(6,6,6);

dRdX0 = dR0_d0(R0);
for i = 1:6
    P0_d0(:,:,i)=[dRdX0(:,:,i) zeros(3); zeros(3) zeros(3) ];
end

