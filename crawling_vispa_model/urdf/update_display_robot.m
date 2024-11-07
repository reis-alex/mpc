function [display]=update_display_robot(robot,q,R0,r0)
%--- Kinematics ---%
%Kinematics
[~,RL,rJ,rL,~,~]=Kinematics(eye(3),zeros(3,1),q,robot);

for ind=1:robot.n_links_joints
    if ~isempty(robot.links(ind).visual.vertices)
        ind_parent=(robot.links(ind).parent_joint);
        T_joint=[RL(:,:,ind_parent) rJ(:,ind_parent); 0 0 0 1];
        T=T_joint*robot.links(ind).visual.T;
        v=[robot.links(ind).visual.vertices ones(length(robot.links(ind).visual.vertices),1)];
        vert=T*v';
        if isfield(get(display.body(ind)),'Visible')
            display.body(ind).Visible='on';
            display.body(ind).Vertices=vert(1:3,:)';
            display.body(ind).Faces=robot.links(ind).visual.faces;
        else
            if strcmp(robot.links(ind).stif.type,'flexible')
                FaceAlpha=0.4;
            else
                FaceAlpha=0.9;
            end
            display.body(ind)=patch('Faces',robot.links(ind).visual.faces,...
                'Vertices',vert(1:3,:)',...
                'FaceColor',       [0.4 0.8 1.0], ...
                'EdgeColor',       'none',        ...
                'FaceLighting',    'gouraud',     ...
                'AmbientStrength', 0.15, ...
                'FaceAlpha',FaceAlpha);
        end
    else
        if isfield(get(display.body(ind)),'Visible')
        display.body(ind).Visible='off';
        end
    end
    
    display.Joint(ind).XData=rJ(1,ind);
    display.Joint(ind).YData=rJ(2,ind);
    display.Joint(ind).ZData=rJ(3,ind);
    
    display.COM(ind).XData=rL(1,ind);
    display.COM(ind).YData=rL(2,ind);
    display.COM(ind).ZData=rL(3,ind);
end

if ~isempty(robot.base_link.visual.vertices)
    v=[robot.base_link.visual.vertices ones(length(robot.base_link.visual.vertices),1)];
    %vert=robot.base_link.visual.T*v';
    vert=[R0, r0; 0 0 0 1]*v';
    display.body(ind)=patch('Faces',robot.base_link.visual.faces,...
        'Vertices',vert(1:3,:)',...
        'FaceColor',       [0.9 0.4 0.4], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15, ...
        'FaceAlpha',.8);
end

drawnow();