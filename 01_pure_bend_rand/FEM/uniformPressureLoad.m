function f = uniformPressureLoad(P,S,nodes,t)
% Construct a nodal force vector f from a uniform pressure load of
% magnitude P on a surface S. Vector f is of same size as the global force
% vector, and can be added to it directly.
% 
% Inputs:
% P:     2-by-1 array of pressure components in x- and y-directions
% S:     n-by-2 or n-by-3 array of surface elements
% nodes: global list of node coordinates
% t:     element thickness (assumed same for all)
%
% Outputs:
% f:     global nodal force vector contribution
%
% Example
%
% To apply a pressure p to the surface patch below (comprising linear els)
%
%       p____________________p             y
%       |  |  |  |  |  |  |  |             |
%  o----o-------o--o-----o---o---------o   |___x
%  1    2       3  4     5   6         7
%
% P and S would be:
% P = [0 -p]
% S = [2 3;3 4;4 5;5 6]
% 
% f would be:
% f = [0 0 0 f2y 0 f3y 0 f4y 0 f5y 0 f6y 0 0]'
% where f2y, f3y, etc. are the only non-zero values.
%
% For quadratic elements, with 3 nodes per edge, the procedure is the same,
% but S is now n-by-3. Each row of S contains node indices [ni,nj,nk],
% where ni & nj are at the ends of the corresponding element, and nk is at
% its mid-point. I.e. the surface element looks like:
% o------o------o
% ni     nk     nj
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2016.

nE = size(S,1); % Number of surface elements

nN = size(nodes,1); % Number of nodes
f = zeros(2*nN,1); % Initialise force vector

if size(S,2) == 2 % Linear element version
   for i = 1:nE
      ni = S(i,1); nj = S(i,2); % Node id's for current surface element
      L = norm(nodes(ni,:)-nodes(nj,:)); % Current element length
      fiX = L*t*P(1); % x-force on this element
      f(2*ni-1) = f(2*ni-1) + 0.5*fiX; % Add contributions to global list
      f(2*nj-1) = f(2*nj-1) + 0.5*fiX;
      fiY = L*t*P(2); % y-force on this element
      f(2*ni) = f(2*ni) + 0.5*fiY; % Add contributions to global list
      f(2*nj) = f(2*nj) + 0.5*fiY; % Add contributions to global list
   end
elseif size(S,2) == 3 % Quadratic element version
    for i = 1:nE
       ni = S(i,1); nj = S(i,2); nk = S(i,3); % Node id's for current surface element
       L = norm(nodes(ni,:)-nodes(nj,:)); % Current element length
       fiX = L*t*P(1); % x-force on this element
       f(2*ni-1) = f(2*ni-1) + fiX/6; % Add contributions to global list
       f(2*nj-1) = f(2*nj-1) + fiX/6;
       f(2*nk-1) = f(2*nk-1) + 2*fiX/3;
       fiY = L*t*P(2); % y-force on this element
       f(2*ni) = f(2*ni) + fiY/6; % Add contributions to global list
       f(2*nj) = f(2*nj) + fiY/6;
       f(2*nk) = f(2*nk) + 2*fiY/3;
    end
else
    warning('Invalid surface definition - returning zero force vector')
end
