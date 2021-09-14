% Parameters
k = 30;  h = 60; Tinf = 25; Q1 = 0; Q2 = 1e6; qn = 2e5;

nele = 2;         % no of elements
nnode = 7;        % no of nodes

% Gauss Points and weights for one point gauss quadrature
% Gausspoints for 2x2 line integration
% --------------------------------
ngauss = 2;

xivec = [-0.57735, 0.57735];
xi_w = [1, 1];

etavec = [-0.57735, 0.57735];
eta_w = [1, 1];


% co-ordinates for elements
% -------------------------
coord = [1, 0.0, 0.0;    % 1st column is element no., 2nd column is x-coordinate, 3rd column is y-coordinate
         2, 0.5, 0.0; 
         3, 0.8, 0.0; 
         4, 0.8, 0.3;
         5, 0.5, 0.3;
         6, 0.5, 0.6;
         7, 0.0, 0.6];

connect = [1, 1, 2, 6, 7;     % 1st column is element no., other columns are global node no. of that particular element
           2, 2, 3, 4, 5];
