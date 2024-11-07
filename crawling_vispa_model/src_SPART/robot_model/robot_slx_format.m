function [robot_sim]= robot_slx_format(robot)

robot_sim.n_q=robot.n_q;

robot_sim.n_links_joints=robot.n_links_joints;
robot_sim.n_links_flex= robot.n_links_flex;

robot_sim.joints=robot.joints;
robot_sim.flexs=robot.flexs;
robot_sim.base_link=robot.base_link;
robot_sim.con=robot.con;

for ind=1:robot_sim.n_links_joints
    
    robot_sim.links(ind).id=robot.links(ind).id;
    robot_sim.links(ind).parent_joint=robot.links(ind).parent_joint;
    robot_sim.links(ind).T=robot.links(ind).T;
    robot_sim.links(ind).mass=robot.links(ind).mass;
    
    robot_sim.links(ind).inertia=robot.links(ind).inertia;
    
%     robot_sim.links(ind).stif.nb_mode=robot.links(ind).stif.nb_mode;
%     robot_sim.links(ind).stif.pulse=robot.links(ind).stif.pulse;
%     robot_sim.links(ind).stif.Li=robot.links(ind).stif.Li;
%     robot_sim.links(ind).stif.damp=robot.links(ind).stif.damp;
    
end