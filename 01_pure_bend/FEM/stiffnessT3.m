function Ke = stiffnessT3(X,t,E,nu,conf)
% Construct the 6-by-6 stiffness matrix for an arbitrarily shaped 3-node
% triangular element. 1-point integration is used.
%
% Inputs:
% - X:    3-by-2 matrix of element nodal coords
% - t:    element thickness
% - E:    Young's modulus
% - nu:   Poisson's ratio
% - conf: plane strain or plane stress configuration
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2016.

% Check inputs
if nargin ~= 5
   error('5 inputs required')
end
[m,n] = size(X);
if m ~= 3 || n ~= 2
   error('X must be 3-by-2')
end

% Quad points and weights (1-pt, natural coords)
Q = [1/3 1/3];
% !! Point coords not actually required in this case, since shape function
% derivs are constant.
W = 0.5;

D = [];
if conf == 1 % Plane strain case
   D = (E/(1+nu)/(1-2*nu)) * [1-nu nu 0;
                              nu 1-nu 0;
                              0 0 0.5*(1-2*nu)];
elseif conf == 2 % Plane stress case
   D = (E/(1-nu*nu)) * [1 nu 0 ;
                     nu 1 0 ;
                     0 0 0.5*(1-nu)];
else
   error('Unknown configuration type')
end

% Stiffness matrix
dNdL = [1  0; % Shape function derivs (nat. coords)
        0  1;
       -1 -1];
J = dNdL'*X; % Jacobian
dNdx = dNdL*inv(J)'; % Shape function derivs (global coords)
B = zeros(3,3*2); % Strain-displacement matrix
for j = 1:3
   B(:,2*j-1:2*j) = [dNdx(j,1),0;
                     0,dNdx(j,2);
                     dNdx(j,2),dNdx(j,1)];
end
Ke = t*B'*D*B*det(J)*W; % Numerical integration to compute k (only 1 point)
