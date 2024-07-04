function [En,Sn] = nodalStrainsStresses(nodes,els,U,E,nu,conf)
% Compute the nodal strains and stresses in a 2D model.
%
% Inputs:
% - nodes: matrix of node coordinates
% - els:   matrix of element connectivities
% - U:     vector of nodal displacements
% - E:     Young's modulus (assumed same for all elements)
% - nu:    Poisson's ratio (assumed same for all elements)
% - conf: plane strain or plane stress configuration
%
% Returns:
% - En:    matrix of nodal strains (n-by-3):
%          row i contains [Exx,Eyy,Gxy] for node i
% - Sn:    matrix of nodal stresses (n-by-3):
%          row i contains [Sxx,Syy,Txy] for node i
%
% Strains and stresses are calculated within elements, but in general, the
% results are not continuous across element boundaries. To produce smooth
% distributions, the values are usually averaged at the nodes.
%
% Essentially, a two step process:
% 1. For each element, compute the strains/stresses at its nodes
% 2. For each node, average the contributions from each element
%
% Currently implemented to support 3-node triangle (T3) elements only.
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2017.

nnds = size(nodes,1);
[nels,npe] = size(els);

% Stress points (evaluation points - i.e. corners) in natural coords:
if npe == 3 % triangle elements
   stressPts = [0 0;1 0;0 1];
elseif npe == 4 % triangle elements  此处修改
   stressPts = [-1/sqrt(3) -1/sqrt(3);1/sqrt(3) -1/sqrt(3);1/sqrt(3) 1/sqrt(3);-1/sqrt(3) 1/sqrt(3)];
else
end

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

% Number of elements connected to each node. We need this to compute the
% average of the contributions from each element.
nodeValence = zeros(nnds,1);

En = zeros(nnds,3);
Sn = zeros(nnds,3);
% Loop over elements, computing strains/stresses at each of their nodes
for i = 1:nels
   % Current element connectivities, DOFs...:
   conn = els(i,:);
   npe = length(conn); % nodes per element
   dofs = zeros(2*npe,1);
   dofs(1:2:end-1) = conn*2-1; % DOF mapping
   dofs(2:2:end) = conn*2;
   % Increment the nodeValence tallies:
   nodeValence(conn) = nodeValence(conn) + 1;
   
   % Current element node coords
   X = nodes(conn,:);
   
   % Loop over stress evaluation points in the element (i.e. corners)
   for q = 1:npe
      r = stressPts(q,1); s = stressPts(q,2); % Natural coords of current point
      % Actually, r & s are redundant for T3 element case, since
      % strain/stress are constant 此处修改
      if npe==3
      dNdr = [-1 -1; % Shape function derivs (nat. coords)
               1  0;
               0  1];        
      J = dNdr'*X; % Jacobian
      dNdx = dNdr*inv(J)'; % Shape function derivs (global coords)
      B = zeros(3,npe*2); % Strain-displacement matrix
      for j = 1:npe % Assemble it in blocks of 2 columns
         B(:,2*j-1:2*j) = [dNdx(j,1),0;
                           0,dNdx(j,2);
                           dNdx(j,2),dNdx(j,1)];
      end
      
      
      elseif npe==4
            syms ksi eta
dNdL = 1/4*[-(1-eta) -(1-ksi); % Shape function derivs (nat. coords)
             (1-eta) -(1+ksi);
             (1+eta)  (1+ksi);
            -(1+eta)  (1-ksi)];
J = dNdL'*X; % Jacobian
dNdx = dNdL*inv(J)'; % Shape function derivs (global coords)
%dNdx=simplify(dNdx);
B =vpa( zeros(3,4*2)); % Strain-displacement matrix
for j = 1:npe
   B(:,2*j-1:2*j) =[dNdx(j,1),0;
                     0,dNdx(j,2);
                     dNdx(j,2),dNdx(j,1)];
end
       ksi=r;eta=s;
       B=eval(B);
      end
      % Strains and stresses at current point
      strain = B*U(dofs);
      stress = D*strain;
      % Add to global tallies
      En(conn(q),:) = En(conn(q),:) + strain';
      Sn(conn(q),:) = Sn(conn(q),:) + stress';
   end
end

% Convert nodal values from sums into averages:
En = En./[nodeValence nodeValence nodeValence];
Sn = Sn./[nodeValence nodeValence nodeValence];