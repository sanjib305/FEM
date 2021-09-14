function [ke,fe] = domain(xi, eta, coord, k, Q)
% This function calculate element Stiffness Matrix

% INPUT:
% ======
% k = Thermal_conductivity
% Q = Heat generation per unit volume
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix
% fe = Element load vector

ddmat = 0.25* [-(1-eta), (1-eta), (1+eta), -(1+eta);
               -(1-xi), -(1+xi),  (1+xi),   (1-xi)];

J = ddmat*coord;

N = 0.25* [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];

B = inv(J)*ddmat;


ke = k*det(J)*(B'*B);

fe = N'*Q*det(J);

end
