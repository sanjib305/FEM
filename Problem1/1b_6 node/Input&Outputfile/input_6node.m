% Parameters
k = 30;  h = 60; Tinf = 25; Q1 = 0; Q2 = 1e6; qn = 2e5;

nele = 4;         % no of elements
nnode = 16;        % no of nodes

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
         2, 0.5, 0.0; 
         3, 0.8, 0.0; 
         4, 0.8, 0.3;
         5, 0.5, 0.3;
         6, 0.5, 0.6;
         7, 0.0, 0.6; 
         8, 0.25, 0.0;
         9, 0.65, 0.0;
        10, 0.8, 0.15;
        11, 0.65, 0.3;
        12, 0.25, 0.6;
        13, 0.0, 0.3;
        14, 0.25, 0.3;
        15, 0.5, 0.15;
        16, 0.65, 0.15];

connect = [1, 1, 2, 7, 8, 14, 13;     % 1st column is element no., other columns are global node no.s of that particular element
           2, 2, 6, 7, 5, 12, 14;
           3, 2, 4, 5, 16, 11, 15;
           4, 2, 3, 4, 9, 10, 16];
