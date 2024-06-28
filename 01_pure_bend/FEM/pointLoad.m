function f = pointLoad(plInd,plVal,nN)
% Construct a nodal force vector corresponding to point loads. Vector f is
% of same size as the global force vector, and can be added to it directly.
% 
% Inputs:
% plInd: n-vector of indices of loaded DOFs
% plVal:	n-vector of corresponding force magnitudes
% nN:    number of nodes in the mesh
%
% Outputs:
% f:     global nodal force vector contribution
%
% Example
%
% To apply the following loads
%
%           P1        P2
%           |         |
%           V         V
% o----o----o----o----o----o----o <--P3
% 1    2    3    4    5    6    7
%
% plInd and plVal would be:
% plInd = [6;10;13]
% plVal = [P1;P2;P3]
%
% I.e. the y-dofs of nodes 3 & 5, and the x-dof of node 7 are loaded with
% forces P1, P2, and P3, respectively.
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2017.

f = zeros(2*nN,1); % Initialise force vector

f(plInd) = plVal; % Set values
