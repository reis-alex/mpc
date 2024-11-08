function [voro] = voronoicell(robot,maxrows,Xc)

[~,~,n_robot] = size(robot);

for i = 1
interp_x{1} = robot(1,:,i);
interp_y{1} = robot(2,:,i);

others = find([1:n_robot]~=i);

for j = 1:length(others)
    interp_x{j+1} = robot(1,:,others(j));
    interp_y{j+1} = robot(2,:,others(j));
    sym{j} = reflect([robot(1,:,i); robot(2,:,i)],[robot(1,:,others(j)); robot(2,:,others(j))]);
end

interp = vertcat(horzcat([interp_x{:}]',[interp_y{:}]'),[sym{:}]');
dt      = delaunayTriangulation(interp);
[V,R_v] = voronoiDiagram(dt);
jj = 1;

% generate box
Vbox = [];
Vbox = [Vbox; [robot(1,:,i)  robot(2,:,i)] + 1*[+0 +0.5]];
Vbox = [Vbox; [robot(1,:,i)  robot(2,:,i)] + 1*[-0 -0.5]];
Vbox = [Vbox; [robot(1,:,i)  robot(2,:,i)] + 1*[+0.5 -0]];
Vbox = [Vbox; [robot(1,:,i)  robot(2,:,i)] + 1*[-0.5 -0]];
box = Polyhedron('V',Vbox);

vert = V([R_v{jj}(:)],:);
voro{i}(jj) = Polyhedron('V',vert(~isinf(vert(:,1)),:));
% voro{i}(jj) = voro{i}(jj).intersect(Xc);
% voro{i}(jj) = voro{i}(jj).intersect(box);
voro{i}(jj).minHRep();

clear interp_x interp_y interp V R_v dt sym vert j jj Vbox box

end
