function [branch,child,child_base]=ConnectivityMap(robot)
% Produces the connectivity map for a robot model.
%
% [branch,child,child_base]=ConnectivityMap(robot)
%
% :parameters:
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * branch -- Branch connectivity map. This is a [nxn] lower triangular matrix. If the i,j element is 1 it means that the ith and jth link are on the same branch. 
%   * child -- A [nxn] matrix. If the i,j element is 1, then the ith link is a child of the jth link.
%   * child_base -- A [nx1] matrix. If the ith element is 1, the ith link is connected to the base-link.
%
% See also: :func:`src.robot_model.urdf2robot` and :func:`src.robot_model.DH_Serial2robot`.

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

%Pre-allocate branch connectivity map
branch=zeros(robot.n_links_joints,robot.n_links_joints);

%Populate branch connectivity map
for i=robot.n_links_joints:-1:1
    last_parent_link=i;
    branch(i,i)=1;
    %Descent through the branch until the base, populating the connectivity map
    while true
        parent_link=robot.joints(robot.links(last_parent_link).parent_joint).parent_link;
        if parent_link==0
            break
        end
        branch(i,parent_link)=1;
        last_parent_link = parent_link;
    end
end

%Populate child map
child=zeros(robot.n_links_joints,robot.n_links_joints);
child_base=zeros(robot.n_links_joints,1);
for i=robot.n_links_joints:-1:1
    parent_link=robot.joints(robot.links(i).parent_joint).parent_link;
    if parent_link~=0
        child(i,parent_link)=1;
    else
        child_base(i)=1;
    end
end

end
