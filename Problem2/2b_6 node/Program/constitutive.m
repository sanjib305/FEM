function [C] = constitutive(E, nu)

% This function calculate Constitutive Matrix

% INPUT:
% ======
% E = Young's Modulus
% nu = Poisson's ratio
%
% OUTPUT:
% =======
% C = Constitutive Matrix

   
C = [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
C = C*E/(1-nu^2);

    
end
