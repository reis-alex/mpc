function DCM = Angles321_DCM(Angles)
% Convert the Euler angles (321 sequence), x-phi, y-theta, z-psi to its DCM equivalent.
%
% DCM = Angles321_DCM(Angles)
%
% :parameters: 
%   * Angles -- Euler angles [x-phi, y-theta, z-psi] -- [3x1].
%
% :return: 
%   * DCM -- Direction Cosine Matrix -- [3x3].
%
% See also: :func:`src.attitude_transformations.Angles123_DCM` and :func:`src.attitude_transformations.DCM_Angles321`.

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

phi = Angles(1);
theta = Angles(2);
psi = Angles(3);

DCM=[   cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
        sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta);
        cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)];


end