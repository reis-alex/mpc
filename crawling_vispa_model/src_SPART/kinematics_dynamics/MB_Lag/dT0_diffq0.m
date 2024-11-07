function dT0_dq0 = dT0_diffq0(q0,R0)

    phi = q0(1);
    theta = q0(2);
    psi = q0(3);
   
    c0z = cos(psi);
    s0z = sin(psi);
    c0y = cos(theta);
    s0y = sin(theta);
    c0x = cos(phi);
    s0x = sin(phi);
        
    dR_rX = zeros(3);
    dR_rY = zeros(3);
    dR_rZ = zeros(3);
   
    %Create individual DCM and derivatives
    C1=[1 0         0;
        0 cos(phi) -sin(phi);
        0 sin(phi) cos(phi)]; 
    
    C2=[cos(theta) 0 sin(theta);
        0          1 0;
       -sin(theta) 0 cos(theta)]; 
    
    C3=[cos(psi) -sin(psi) 0;
        sin(psi)  cos(psi) 0;
        0         0        1];
    
    

%     dC1 = SkewSym([1;0;0])*C1;
%     dC2 = SkewSym([0;1;0])*C2;
%     dC3 = SkewSym([0;0;1])*C3;


    dC1=[0  0         0;
            0 -sin(phi) -cos(phi);
            0  cos(phi) -sin(phi)]; 
    
    dC2=[-sin(theta) 0  cos(theta);
               0          0  0;
              -cos(theta) 0 -sin(theta)]; 
    
    dC3=[-sin(psi) -cos(psi) 0; 
             cos(psi) -sin(psi) 0;
             0         0        0];

    %Compute DCM derivatives for 123 sequence
    
    DCM_phi   = C3*C2*dC1;
    DCM_theta = C3*dC2*C1;
    DCM_psi   = dC3*C2*C1;

    dRdq0 = zeros(3,3,3, "like", R0);
    dRdq0(:,:,1) = DCM_phi ;
    dRdq0(:,:,2) = DCM_theta;
    dRdq0(:,:,3) = DCM_psi;


     dT0_dq0 = zeros(4,4,6, "like", R0);
%     dT0dq0(:,:,1) = [dRdq0(:,:,1) zeros(3,1);zeros(1,4)];
%     dT0dq0(:,:,2) = [dRdq0(:,:,2)  zeros(3,1);zeros(1,4)];
%     dT0dq0(:,:,3) = [dRdq0(:,:,3) zeros(3,1) ;zeros(1,4)]; 
%     dT0dq0(:,:,4) = [dRdq0(:,:,4), [1 0 0]';zeros(1,4)];
%     dT0dq0(:,:,5) = [dRdq0(:,:,5), [0 1 0]';zeros(1,4)];  
%     dT0dq0(:,:,6) = [dRdq0(:,:,6), [0 0 1]';zeros(1,4)]; 

for i = 1:3
    dT0_dq0(1:3,1:3,i) = dRdq0(:,:,i);%SkewSym(R0(:,i))'*R0;
    dT0_dq0(i,4,i+3)   = 1; 

end




end