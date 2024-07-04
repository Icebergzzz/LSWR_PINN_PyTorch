function K = stiffnessGlobal(nodes,els,t,E,nu,conf)
% Assemble the global stiffness matrix for a 2D FE problem.
%
% Inputs:
% - nodes: global list of node coordinates
% - els:   element connectivity array
% - t:     thickness of elements
% - E:     Young's modulus
% - nu:    Poisson's ratio
% - conf:  plane strain or plane stress configuration
%
% NB: single values of t, E and nu should be provided, i.e. these are
% assumed constant for all elements in the model.
% 
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2017.

[nels,npe] = size(els);
%if (npe ~= 3)
%   error('Only T3 elements supported, so far')
%end
nnds = size(nodes,1);

% Initialise K using sparse representation - much more memory efficient,
% and essential for larger models:
K = sparse(2*nnds,2*nnds);

% Loop on elements, construct element matrices, assemble into global
for i = 1:nels
   conn = els(i,:); % Connectivity of current element
   X = nodes(conn,:); % Array of element node coords
   if npe == 3
      Ke = stiffnessT3(X,t,E,nu,conf);
   elseif npe == 4 
      Ke = stiffQ4(X,t,E,nu,conf); % Place holder for Q4 element code
      else
 end
   dofs = zeros(2*npe,1);
   dofs(1:2:end-1) = conn*2-1; % DOF mapping
   dofs(2:2:end) = conn*2;
   K(dofs,dofs) = K(dofs,dofs) + Ke; % Add to global stiffness matrix
end
