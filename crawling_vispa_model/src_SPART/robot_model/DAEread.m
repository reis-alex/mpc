function [F,V]=DAEread(filename)
% DAEread imports geometry from an DAE file into MATLAB.
%
%    [F,V] = DAEREAD(FILENAME) returns the faces F and vertices V separately.

%Read DAE file
DAE=xmlread(filename);

%Only one collada object description per DAE file is allowed
collada_DAE = DAE.getElementsByTagName('COLLADA');
if collada_DAE.getLength() ~= 1
    error('URDF contains more than one collada object. This function only accepts a single collada object descriptioon per xml');
end
collada_DAE =collada_DAE.item(0);


%Get geometry
geometries = collada_DAE.getElementsByTagName('library_geometries');
geometry = geometries.item(0).getElementsByTagName('geometry');
%geometry_id= char(geometry.item(0).getAttribute('id'));
mesh= geometry.item(0).getElementsByTagName('mesh');

position=mesh.item(0).getElementsByTagName('source').item(0);
%normal=mesh.item(0).getElementsByTagName('source').item(1);
triangles=mesh.item(0).getElementsByTagName('triangles').item(0);

fa_pos=position.getElementsByTagName('float_array').item(0);
pos_array=eval([ '[ ' char(fa_pos.getFirstChild.getData) ']' ]);
nb_pos= eval(position.getElementsByTagName('accessor').item(0).getAttribute('count'));

V=zeros(nb_pos,3);
V(:,1:3)=(reshape(pos_array,3,nb_pos))';
%fa_norm=normal.getElementsByTagName('float_array').item(0);
%norm_array=eval([ '[ ' char(fa_norm.getFirstChild.getData) ']' ]);

ver_tri=triangles.getElementsByTagName('p').item(0);
ver_array=eval([ '[ ' char(ver_tri.getFirstChild.getData) ']' ]);

nb_face=eval(char(triangles.getAttribute('count')));

face_norm=reshape(ver_array,2,3*nb_face);


face=face_norm(1,:);

F=zeros(nb_face,4);
F(:,1:3)=(reshape(face,3,nb_face))';
F(:,4)=F(:,1);
F=F+ones(nb_face,4);





