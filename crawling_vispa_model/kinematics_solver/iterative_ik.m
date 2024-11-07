function [q,e]=iterative_ik(Rd,rd,qm,robot,ind,range)
display=0;
step_size=0.01;
tolerance=0.001;
lambda=1e-12;
if display
    display_handle=display_robot(robot,eye(3),[0;0;0],qm);
    pause
end

[Rp,rp]=mgd_rel(qm,ind,robot);
e = delta_pose(Rd,rd,Rp,rp);
q=qm;
it =0;
prev_e=e;
while norm(e) >= tolerance && it<1000  %&& norm(prev_e(4:6))>=norm(e(4:6))
    prev_e=e;
    [Jp,~,~]=jacob_rel(q,ind,robot,range);
    
    %e(1:3,1)=e(1:3,1)*0.01;
    delta_q = (Jp' * Jp + lambda * eye(6))\(Jp'*e);
    %if norm(delta_q)>step_size
        q(range,1) = q(range,1)+ step_size * delta_q;
    % else
    %     q(range,1) = q(range,1)+delta_q;
    % end
    [Rp,rp]=mgd_rel(q,ind,robot);
    e = delta_pose(Rd,rd,Rp,rp);
    if display
        update_display_robot(display_handle,robot,eye(3),[0;0;0],q);
        drawnow();
    end
    fprintf('att= %f\t pos= %f,%f, delta_q= %f cond= %f\n',norm(e(1:3)),norm(e(4:6)),norm(prev_e(4:6)),norm(delta_q),cond(Jp));
    disp(e(1:3)');
    it=it+1;
    %pause
end
end