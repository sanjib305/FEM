function [fe] = gamat(beta, coord, P0, P1, theta, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% tvec = traction components
% xi = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge traction is prescribed
%
% OUTPUT:
% =======
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi = (1+beta)/2;

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N1 = 1-xi; N2 = xi; N3 = 0;
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N1 = 0; N2 = 1-xi; N3 = xi;
elseif edgeno == 3
    ledge = sqrt((x(4)-x(3))^2 + (y(4)-y(3))^2);
    N1 = xi; N2=0; N3 = 1-xi;
end

N = zeros(2,6);
N(1,1:2:6) = [N1, N2, N3];
N(2,2:2:6) = [N1, N2, N3];

theta = theta*(pi/180);
y = N1*y(1) + N2*y(2) + N3*y(3);
tvec= [P1*cos(theta) + (P0 - P1)*cos(theta)*y/0.4;
       P1*sin(theta) + (P0 - P1)*sin(theta)*y/0.4];

fe = ledge*N'*tvec*0.5;
  
end
