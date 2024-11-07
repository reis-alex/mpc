function dRdX0 = dR0_d0(R0) %[dR0_X, dR0_Y, dR0_Z, dR_rX, dR_rY, dR_rZ] = dR0_d0(theta0)
%DR0_D0 Return derivatives of euler XYZ rotation respect to base_link rotation theta0 and translation (always zero matrix when derivating over translation) 
%   Detailed explanation goes here
    
    theta0 = DCM_Angles321(R0);

    theta_x = theta0(1);
    theta_y = theta0(2);
    theta_z = theta0(3);
   
    c0z = cos(theta_z);
    s0z = sin(theta_z);
    c0y = cos(theta_y);
    s0y = sin(theta_y);
    c0x = cos(theta_x);
    s0x = sin(theta_x);
    
    
    drX12 =  c0z*s0y*c0x + s0z*s0x;
    drX13 = -c0z*s0y*s0x;
    drX22 =  c0z*c0x + s0z*s0y*c0x;
    drX23 = -c0z*c0x -s0z*s0y*s0x;
    drX32 =  c0y*c0x;
    drX33 = -c0y*s0x;

    dR0_X = [0     drX12 drX13;
             0     drX22 drX23;
             0     drX32 drX33];

    drY11 = -c0z*s0y;
    drY12 =  c0z*c0y*s0x;
    drY13 =  s0z*c0y + c0z*c0y*c0x;
    drY21 = -s0z*s0y;
    drY22 =  s0z*c0y*s0x;
    drY23 =  s0z*c0y*c0x;
    drY31 = -c0y;
    drY32 = -s0y*s0x;
    drY33 = -s0y*c0x;

    dR0_Y = [drY11 drY12 drY13;
             drY21 drY22 drY23;
             drY31 drY32 drY33];

    drZ11 = -s0z*c0y;
    drZ12 = -s0z*s0y*s0x - c0z*c0x;
    drZ13 =  c0z*s0y - s0z*s0y*s0x;
    drZ21 =  c0z*c0y;
    drZ22 = -s0z*s0x + c0z*s0y*s0x;
    drZ23 =  s0z*s0x + c0z*s0y*c0x;

    dR0_Z = [drZ11 drZ12 drZ13;
             drZ21 drZ22 drZ23;
             0     0     0   ];
    
         
    dR_rX = zeros(3);
    dR_rY = zeros(3);
    dR_rZ = zeros(3);
    
    dRdX0 = zeros(3,3,6, "like", R0);
    dRdX0(:,:,1) = dR0_X;
    dRdX0(:,:,2) = dR0_Y;
    dRdX0(:,:,3) = dR0_Z;
    dRdX0(:,:,4) = dR_rX;
    dRdX0(:,:,5) = dR_rY;
    dRdX0(:,:,6) = dR_rZ;
end

