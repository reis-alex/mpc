function [robot,robot_keys]=urdf2robot(filename,verbose_flag)
% Creates a SPART robot model from a URDF file.
%
% [robot,robot_keys] = urdf2robot(filename,verbose_flag)
% 
% :parameters: 
%   * filename -- Path to the URDF file.
%   * verbose_flag -- True for verbose output (default False).
%
% :return:
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%   * robot_keys -- Links/Joints name map (see :doc:`/Tutorial_Robot`).
%
% This function was inspired by:
% https://github.com/jhu-lcsr/matlab_urdf/blob/master/load_ne_id.m

%{  
    LICENSE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

%=== CODE ===%

%Verbose flag defaults to false
if nargin<2
    verbose_flag=false;
end

%Create Robot structure.
robot=struct();

%Read URDF file
urdf=xmlread(filename);

%Only one robot description per URDF file is allowed
robot_urdf = urdf.getElementsByTagName('robot');
if robot_urdf.getLength() ~= 1
    error('URDF contains more than one robot. This function only accepts a single robot descriptioon per URDF');
end
robot_urdf =robot_urdf.item(0);

%Get robot name
robot.name = char(robot_urdf.getAttribute('name'));

%Get links and joints from URDF file
links_urdf = robot_urdf.getElementsByTagName('link');
joints_urdf = robot_urdf.getElementsByTagName('joint');

%Find and remove links elements that are not directly under the robot element
i=0;
while i<links_urdf.getLength
    if ~strcmp(links_urdf.item(i).getParentNode.getNodeName(),'robot')
        links_urdf.item(i).getParentNode.removeChild(links_urdf.item(i));
    else
        i=i+1;
    end
end

%Find and remove joint elements that are not directly under the robot element
i=0;
while i<joints_urdf.getLength
    if ~strcmp(joints_urdf.item(i).getParentNode.getNodeName(),'robot')
        joints_urdf.item(i).getParentNode.removeChild(joints_urdf.item(i));
    else
        i=i+1;
    end
end

%Create number of joint variables (to be poulated later)
robot.n_q=[];

%Count links and joints
robot.n_links_joints = links_urdf.getLength;

%Create temporary link and joint maps
links = containers.Map();
joints = containers.Map();

%Display data
if verbose_flag; fprintf('Number of links: %d (including the base link)\n', robot.n_links_joints); end

%Iterate over links
for k = 0:robot.n_links_joints-1
    %Create basic structure with default values
    link = struct();
    link_xml = links_urdf.item(k);
    link.name = char(link_xml.getAttribute('name'));
    link.T=[eye(3),zeros(3,1);zeros(1,3),1];
    link.parent_joint = {};
    link.child_joint = {};
    
    %Grab inertial properties
    inertial = link_xml.getElementsByTagName('inertial').item(0);
    
    %Grab origin properties
    origin = inertial.getElementsByTagName('origin').item(0);
    if ~isempty(origin)
        if ~isempty(char(origin.getAttribute('xyz')))
            link.T(1:3,4) = eval(['[',char(origin.getAttribute('xyz')),']'])';
        end
        if ~isempty(char(origin.getAttribute('rpy')))
            rpy = eval(['[',char(origin.getAttribute('rpy')),']']);
            link.T(1:3,1:3)=Angles321_DCM(rpy')';
        end
    end
    
    %Mass
    mass = inertial.getElementsByTagName('mass').item(0);
    link.mass = eval(char(mass.getAttribute('value')));
    
    %Inertia
    inertia = inertial.getElementsByTagName('inertia').item(0);
    ixx = eval(inertia.getAttribute('ixx'));
    iyy = eval(inertia.getAttribute('iyy'));
    izz = eval(inertia.getAttribute('izz'));
    ixy = eval(inertia.getAttribute('ixy'));
    iyz = eval(inertia.getAttribute('iyz'));
    ixz = eval(inertia.getAttribute('ixz'));
    link.inertia = [ixx, ixy, ixz; ixy, iyy, iyz; ixz, iyz, izz];
    
    %Store this link in the links map
    links(char(link.name))=link;
end

%Iterate over joints
for k = 0:robot.n_links_joints-2  
    %Create basic structure with default values
    joint = struct();
    joint_xml = joints_urdf.item(k);
    joint.name = char(joint_xml.getAttribute('name'));
    joint.type_name = char(joint_xml.getAttribute('type'));
    joint.parent_link = '';
    joint.child_link = '';
    joint.T=[eye(3),zeros(3,1);zeros(1,3),1];
    
    if strcmp(joint.type_name,'revolute') || strcmp(joint.type_name,'continuous')
        joint.type=1;
    elseif strcmp(joint.type_name,'prismatic')
        joint.type=2;
    elseif strcmp(joint.type_name,'fixed')
        joint.type=0;
        joint.axis = [0; 0; 0];
    else
        error(sprintf('Joint type ''%s'' not supported.',joint.type_name));
    end
    
    
    %Get origin properties
    origin = joint_xml.getElementsByTagName('origin').item(0);
    if ~isempty(origin)
        if ~isempty(char(origin.getAttribute('xyz')))
            joint.T(1:3,4) = eval(['[',char(origin.getAttribute('xyz')),']'])';
        end
        if ~isempty(char(origin.getAttribute('rpy')))
            rpy = eval(['[',char(origin.getAttribute('rpy')),']']);
            joint.T(1:3,1:3)=Angles321_DCM(rpy')';
        end
    end
    
    %Get rotation/sliding axis
    axis = joint_xml.getElementsByTagName('axis').item(0);
    if ~isempty(axis)
        joint.axis = eval(['[',char(axis.getAttribute('xyz')),']'])';
    elseif isempty(axis) && joint.type~=0
        %Moving joints need a rotation/sliding axis.
        error_message=[joint.name,' is a moving joint and requires a joint axis.'];
        error(error_message);
    end
    
    %Get parent link name
    parent = joint_xml.getElementsByTagName('parent').item(0);
    if ~isempty(parent)
        joint.parent_link = char(parent.getAttribute('link'));
        
        %Store the joint name in the parent link
        parent=links(joint.parent_link);
        parent.child_joint(end+1) = {joint.name};
        links(joint.parent_link) = parent;
    end
    
    %Get child link name
    child = joint_xml.getElementsByTagName('child').item(0);
    if ~isempty(child)
        joint.child_link = char(child.getAttribute('link'));
        
        %Store the joint name in the child link
        child =links(joint.child_link);
        child.parent_joint(end+1) = {joint.name};
        links(joint.child_link) = child;
    end
    
    %Correct homogeneous transformation so that it is from previous link
    %inertial
    joint.T=parent.T\joint.T;
    
    %Store this joint in the joints map
    joints(char(joint.name))=joint;
end

% Find the base link
for link_name = links.keys
    if isempty(links(char(link_name)).parent_joint)
        base_link = char(link_name);
        if verbose_flag; fprintf('Base link: %s\n',base_link); end
    end
end
%There needs to be a root link
if ~exist('base_link','var')
    error('Robot has no single base link!');
end

%Structure links and joints map into a structure and create a map with
%names and IDs.

%Create ID maps
robot_keys.link_id=containers.Map();
robot_keys.joint_id=containers.Map();
robot_keys.q_id=containers.Map();

%Remove base link from the number of total links
robot.n_links_joints=robot.n_links_joints-1;

%Create links and joints stucture
if robot.n_links_joints>0
    robot.links(robot.n_links_joints) = struct();
    robot.joints(robot.n_links_joints) = struct();
end

%Save base link on its own structure
clink=links(base_link);
robot.base_link.mass=clink.mass;
robot.base_link.inertia=clink.inertia;

%Assign base ID
robot_keys.link_id(base_link)=0;

%Add links and joints into the structure with the standard numbering
nl=-1; %Link index
nj=-1; %Joint index
nq=1; %Joint variable index 
%Recursively scan through the tree structure
for n=1:length(clink.child_joint)
    [robot,robot_keys,nl,nj,nq]=urdf2robot_recursive(robot,robot_keys,links,joints,joints(clink.child_joint{n}),nl+1,nj+1,nq);
end

%Populate number of joint variables
robot.n_q=nq-1;
if verbose_flag; fprintf('Number of joint variables: %d\n',robot.n_q); end;

%--- Add Conectivity Map ---%
[branch,child,child_base]=ConnectivityMap(robot);
robot.con.branch=branch;
robot.con.child=child;
robot.con.child_base=child_base;

end

%--- Recursive function ---%
function [robot,robot_keys,nl,nj,nq]=urdf2robot_recursive(robot,robot_keys,links,joints,child_joint,nl,nj,nq)%#codegen

%Copy the elements of child joint
robot.joints(nj+1).id=nj+1;
robot.joints(nj+1).type=child_joint.type;
%Assign joint variable if joint is revolute or prismatic
if child_joint.type
    robot.joints(nj+1).q_id=nq;
    robot_keys.q_id(child_joint.name)=nq;
    nq=nq+1;
else
    %Fixed joint assign -1
    robot.joints(nj+1).q_id=-1;
end
robot.joints(nj+1).parent_link=robot_keys.link_id(child_joint.parent_link);
robot.joints(nj+1).child_link=nl+1;
robot.joints(nj+1).axis=child_joint.axis;
robot.joints(nj+1).T=child_joint.T;

%Copy elements of child link
clink=links(child_joint.child_link);
robot.links(nl+1).id=nl+1;
robot.links(nl+1).parent_joint=nj+1;
robot.links(nl+1).T=clink.T;
robot.links(nl+1).mass=clink.mass;
robot.links(nl+1).inertia=clink.inertia;

%Assign ID
robot_keys.joint_id(child_joint.name)=nj+1;
robot_keys.link_id(clink.name)=nl+1;

%Recursively scan through the tree structure
for n=1:length(clink.child_joint)
    [robot,robot_keys,nl,nj,nq]=urdf2robot_recursive(robot,robot_keys,links,joints,joints(clink.child_joint{n}),nl+1,nj+1,nq);
end

end




