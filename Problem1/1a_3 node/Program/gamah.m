function [ke,fe] = gamah(beta, coord, h, Tinf, edgeno)
% This function calculate element stiffness matrix and loadvector on gamah part

% INPUT:
% ======
% h = Convective heat transfer coefficient
% Tinf = Ambient temperature
% beta = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge convective heat transfer is taking place
%
% OUTPUT:
% =======
% ke = Element stiffness matrix
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi = (1+beta)/2;

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N = [1-xi, xi, 0];
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N = [0, 1-xi, xi];
elseif edgeno == 3
    ledge = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
    N = [xi, 0, 1-xi];
end

fe = N'*h*ledge*Tinf*0.5;
ke = h*(N'*N)*ledge*0.5;

end