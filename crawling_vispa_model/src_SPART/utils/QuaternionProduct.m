function [q] = QuaternionProduct(q1,q2)
% Computes the product between two quaternions q=q1?q2
%
% Quaternion convention q=[q_vector;q_scalar] (4x1 column vector)
% --- SCALAR PART LAST ---
%
% Inputs:
%   q1 -> quaternion
%   q2 -> quaternion
%
% Outputs:
%   q -> quaternion product q=q1?q2

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

q=Matrix_QuatProd_Left(q1)*q2;

end
