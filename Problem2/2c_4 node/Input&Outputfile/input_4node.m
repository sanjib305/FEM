% Parameters
E = 200e9; nu = 0.3; P0 = 2e3; P1 = 4e3; theta = 20;
nele = 2;         % no of elements
nnode = 6;        % no of nodes

% Gausspoints for 2x2 line integration
% ------------------------------------
ngauss = 2;

xivec = [-0.57735, 0.57735];
xi_w = [1, 1];

etavec = [-0.57735, 0.57735];
eta_w = [1, 1];


% co-ordinates for elements
% -------------------------
coord = [1, 0.0, 0.0;    % 1st column is element no., 2nd column is x-coordinate, 3rd column is y-coordinate
         2, 0.3, 0.0; 
         3, 0.6, 0.0; 
         4, 0.6, 0.4;
         5, 0.3, 0.4;
         6, 0.0, 0.4];

connect = [1, 1, 2, 5, 6;     % 1st column is element no., other columns are global node no.s of that particular element
           2, 2, 3, 4, 5];

