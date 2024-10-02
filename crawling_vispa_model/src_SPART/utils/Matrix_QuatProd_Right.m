function [q_otimes_Mat_Right]=Matrix_QuatProd_Right(q)
% Computes the matrix form of the quaternion product p?q=[q]p
%
% Quaternion convention q=[q_vector;q_scalar] (4x1 column vector)
% --- SCALAR PART LAST ---
%
% Inputs:
%   q -> Quaternion
%
% Outputs:
%   q_otimes_Mat -> Matricial form of the quaternion product p?q=[q]p

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

%=== CODE ===%

q_otimes_Mat_Right=[ q(4), q(3),-q(2), q(1);
                    -q(3), q(4), q(1), q(2);
                     q(2),-q(1), q(4), q(3);
                    -q(1),-q(2),-q(3), q(4)];

end