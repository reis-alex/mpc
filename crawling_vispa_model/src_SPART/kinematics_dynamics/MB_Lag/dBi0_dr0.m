function dBi0_X0 = dBi0_dr0()
%DBI0_DR0 Summary of this function goes here
%   Detailed explanation goes here

dBi0_X0 = zeros(6,6,6);

dBi0_dthetx = zeros(6);
dBi0_dthety = zeros(6);
dBi0_dthetz = zeros(6);

dBi0_dX =  [zeros(3),   zeros(3);
           [0  0  0;
            0  0 -1;
            0  1  0],   zeros(3)];

dBi0_dY =  [zeros(3),   zeros(3);
           [0  0  1;
            0  0  0;
           -1  0  0],   zeros(3)];

dBi0_dZ =  [zeros(3),   zeros(3);
           [0  1  0;
           -1  0  0;
            0  0  0],   zeros(3)];



dBi0_X0(:,:,1) = dBi0_dthetx;
dBi0_X0(:,:,2) = dBi0_dthety;
dBi0_X0(:,:,3) = dBi0_dthetz;


dBi0_X0(:,:,4) = dBi0_dX;
dBi0_X0(:,:,5) = dBi0_dY;
dBi0_X0(:,:,6) = dBi0_dZ;

end

