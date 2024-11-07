function [display]=display_robot(robot,q,R0,r0)
% Display a SPART robot model created from a URDF file.
w=get(0,'ScreenSize');

figure('Position',[0 0 w(3)*2/3 w(4)*2/3]),
% figure()
hold on

if nargin < 4
    R0 = eye(3);
    r0 = zeros(3,1);
    disp("default")
end

%--- Kinematics ---%
%Kinematics
%[RJ,RL,rJ,rL,e,g]=Kinematics(eye(3),zeros(3,1),q,robot);
q0 = [DCM_Angles321(R0);r0];
[dTJ0,dTL0,dTJ,dTL,RJ,RL,rJ,rL,e,g]=Kinematics_and_diff(R0,r0,q0,q,robot);

for ind=1:robot.n_links_joints
    if ~isempty(robot.links(ind).visual.vertices)
        ind_parent=(robot.links(ind).parent_joint);
        T_joint=[RL(:,:,ind_parent) rJ(:,ind_parent); 0 0 0 1];
        T=T_joint*robot.links(ind).visual.T;
        v=[robot.links(ind).visual.vertices ones(length(robot.links(ind).visual.vertices),1)];
        vert=T*v';
%         if strcmp(robot.links(ind).stif.type,'flexible')
%             FaceAlpha=0.4;
%         else
%             FaceAlpha=0.9;
%         end
        FaceAlpha=0.9;
        display.body(ind)=patch('Faces',robot.links(ind).visual.faces,...
            'Vertices',vert(1:3,:)',...
            'FaceColor',       [0.4 0.8 1.0], ...
            'EdgeColor',       'none',        ...
            'FaceLighting',    'gouraud',     ...
            'AmbientStrength', 0.15, ...
            'FaceAlpha',FaceAlpha);
    end
    display.COM(ind)=plot3(rL(1,ind),rL(2,ind),rL(3,ind),'LineWidth',2,'Color','k','Marker','+');
    display.Joint(ind)=plot3(rJ(1,ind),rJ(2,ind),rJ(3,ind),'LineWidth',2,'Color','r','Marker','o');
    if ind==1
        camlight('headlight');
    end
    material('dull');
    axis('image');
    view([-135 35]);
end
ind=ind+1;
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

xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
grid on