function dpose=delta_pose(Rd,rd,Rm,rm)
dang=DCM_deltang(Rm'*Rd); %dang is express in the end-effector base!
dpose(1:3,1)=dang';
dpose(4:6,1)=rd-rm;
end
