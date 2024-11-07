function [q,e]=iterative_ik_pos(rd1,rd2,qm,robot,ind,range)
display=0;
step_size=0.01;
tolerance=0.00001;
lambda=1e-6;
if display
    display_handle=display_robot(robot,eye(3),zeros(3,1),qm);
    pause
end
[~,rp1]=mgd_rel(qm,ind-1,robot);
[~,rp2]=mgd_rel(qm,ind,robot);
e=[rd1-rp1;rd2-rp2];
prev_e=4*e;
q=qm;
it =0;
while norm(e) >= tolerance && it<1000 && norm(prev_e)>norm(e)
    prev_e=e;
    [Jm1,~,~]=jacob_rel(q,ind-1,robot,range);
    [Jm2,~,~]=jacob_rel(q,ind,robot,range);
    Jp=[Jm1(4:6,:);Jm2(4:6,:)];
    %e(1:3,1)=e(1:3,1)*0.01;
    delta_q = (Jp' * Jp + lambda * eye(5))\(Jp'*e);
    %if norm(delta_q)>step_size
        q(range,1) = q(range,1)+ step_size * delta_q;
    % else
    %     q(range,1) = q(range,1)+delta_q;
    % end
    [~,rp1]=mgd_rel(q,ind-1,robot);
    [~,rp2]=mgd_rel(q,ind,robot);
    e=[rd1-rp1;rd2-rp2];
    if display
        update_display_robot(display_handle,robot,eye(3),zeros(3,1),q);
        drawnow();
    end
    %fprintf('pos= %f, cond= %f\n',norm(e),cond(Jp));
    %disp(e');
    it=it+1;
    %pause
end
end