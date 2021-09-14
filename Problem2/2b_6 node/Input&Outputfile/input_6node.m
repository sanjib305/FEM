% Parameters
E = 200e9; nu = 0.3; P0 = 2e3; P1 = 4e3; theta = 20;
nele = 2;         % no of elements
nnode = 9;        % no of nodes

% Gauss Points and weights for one point gauss quadrature
% for triangular area integration
% -------------------------------------------------------
ngaussd = 4;
xivec = [1/3, 0.6, 0.2, 0.2];
etavec = [1/3, 0.2, 0.6, 0.2];
wvec = [-27/48, 25/48, 25/48, 25/48];

% Gausspoints for line integration
% --------------------------------
ngaussb = 3;
betavec = [-0.774597, 0, 0.774597];
ewvec = [5/9, 8/9, 5/9];

% co-ordinates for elements
% -------------------------
coord = [1, 0.0, 0.0;    % 1st column is element no., 2nd column is x-coordinate, 3rd column is y-coordinate
         2, 0.6, 0.0; 
         3, 0.6, 0.4; 
         4, 0.0, 0.4;
         5, 0.3, 0.0;
         6, 0.6, 0.2;
         7, 0.3, 0.4;
         8, 0.0, 0.2;
         9, 0.3, 0.2];

connect = [1, 1, 2, 4, 5, 9, 8;     % 1st column is element no., other columns are global node no.s of that particular element
           2, 2, 3, 4, 6, 7, 9];

