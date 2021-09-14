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

alpha = 1 - xi - eta ;
ddmat = [4*xi - 1, 0, 1 - 4*alpha, 4*eta, -4*eta, 4-8*xi-4*eta;
         0,   4*eta - 1,  1-4*alpha, 4*xi, 4-8*eta-4*xi, -4*xi];

J = ddmat*coord;

R1 = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];
  
R2 = zeros(4,4);  R2(1:2,1:2) = inv(J);  R2(3:4,3:4) = inv(J);

R3 = zeros(4,12);
R3(1,1:2:12) = ddmat(1,1:6);
R3(2,1:2:12) = ddmat(2,1:6);
R3(3,2:2:12) = ddmat(1,1:6);
R3(4,2:2:12) = ddmat(2,1:6);

B = R1*R2*R3;

ke = 0.5*det(J)*(B'*C*B);

end