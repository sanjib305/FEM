function [ke] = elestiff(xi, eta, coord, C)

% This function calculate element Stiffness Matrix

% INPUT:
% ======
% C = Constitutive Matrix
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix

x = coord(:,1);
y = coord(:,2);

J = [x(1)-x(3), y(1)-y(3); x(2)-x(3), y(2)-y(3)];

R1 = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];
  
R2 = zeros(4,4);  R2(1:2,1:2) = inv(J);  R2(3:4,3:4) = inv(J);

R3 = [1, 0, 0, 0, -1, 0;
      0, 0, 1, 0, -1, 0;
      0, 1, 0, 0, 0 ,-1;
      0, 0, 0, 1, 0 ,-1];

B = R1*R2*R3;

ke = 0.5*det(J)*(B'*C*B);
end