function Euler = DCM_Euler(DCM)
% Compute the rotation axis and angle from a DCM.
%
% Input:
%   DCM -> Direction Cosine Matrix.
% Output:
%   Euler = [e;alpha] -> e - Rotation axis, alpha - Rotation angle.

%=== LICENSE ===%

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%=== Code ===%

%Rotation Angle
alpha = real(acos(complex((trace(DCM)-1)/2)));
%Rotation Axis
e=1/(2*sin(alpha))*[DCM(2,3)-DCM(3,2);DCM(3,1)-DCM(1,3);DCM(1,2)-DCM(2,1)];

%Compact axis and angle in a single 4x1 vector
Euler = [e;alpha];

end