function q = DCM_quat(DCM) %#codegen
%Convert Direction Cosine Matrix (DCM) to quaterion q
%
% q = [q1;q2;q3;q4] -> With q4 being the scalar part of the quaternion.


%Compute the squares of the quaternion
qs=[1/4*(1+2*DCM(1,1)-trace(DCM));
    1/4*(1+2*DCM(2,2)-trace(DCM));
    1/4*(1+2*DCM(3,3)-trace(DCM));
    1/4*(1+trace(DCM))];

%Find biggest square value
[~,I]=max(qs);

%Determine quaternions
if I==1
    ql=sqrt(qs(1));
    q=[ ql;
        (DCM(1,2)+DCM(2,1))/(4*ql);
        (DCM(3,1)+DCM(1,3))/(4*ql);
        (DCM(2,3)-DCM(3,2))/(4*ql)];
elseif I==2
    ql=sqrt(qs(2));
    q=[ (DCM(1,2)+DCM(2,1))/(4*ql);
        ql;
        (DCM(2,3)+DCM(3,2))/(4*ql);
        (DCM(3,1)-DCM(1,3))/(4*ql)];
elseif I==3
    ql=sqrt(qs(3));
    q=[ (DCM(3,1)+DCM(1,3))/(4*ql);
        (DCM(2,3)+DCM(3,2))/(4*ql);
        ql;
        (DCM(1,2)-DCM(2,1))/(4*ql)];
elseif I==4
    ql=sqrt(qs(4));
    q=[ (DCM(2,3)-DCM(3,2))/(4*ql);
        (DCM(3,1)-DCM(1,3))/(4*ql);
        (DCM(1,2)-DCM(2,1))/(4*ql);
        ql];
else
    q=zeros(4,1);
end

end