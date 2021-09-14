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

ddmat = 0.25* [-(1-eta), (1-eta), (1+eta), -(1+eta);
               -(1-xi), -(1+xi),  (1+xi),   (1-xi)];

J = ddmat*coord;

R1 = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];
  
R2 = zeros(4,4);  R2(1:2,1:2) = inv(J);  R2(3:4,3:4) = inv(J);

R3 = zeros(4,8);
R3(1,1:2:8) = ddmat(1,1:4);
R3(2,1:2:8) = ddmat(2,1:4);
R3(3,2:2:8) = ddmat(1,1:4);
R3(4,2:2:8) = ddmat(2,1:4);

B = R1*R2*R3;
ke = det(J)*(B'*C*B);

end