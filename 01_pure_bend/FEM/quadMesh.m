function [nodes,els] = quadMesh(m,n,L,D)
% Create a regular quadrilateral mesh of a rectangular region. The
% dimensions of the region are L wide and D high. m and n elements along
% the width and height will be used, respectively.
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, April 2013.

% Create nodes
nodes = zeros((m+1)*(n+1),2);
xlen = L/m; % Element x- and y-lengths
ylen = D/n;
counter = 0; % This will keep track of the node ids as we create them
% Since the mesh forms a grid pattern in this case, we can talk about rows
% and columns of nodes and elements.
for j = 1:n+1 % Step through nodes in each row of the mesh
   for i = 1:m+1 % Step through columns
      counter = counter + 1;
      % Position of the current node, which is in column i, row j:
      nodes(counter,:) = [(i-1)*xlen,(j-1)*ylen - D/2];
   end
end

% Create elements
els = zeros(m*n,4);
counter = 0; % Now, counter keeps track of the element ids
for j = 1:n
   for i = 1:m
      counter = counter + 1;
      % Need, first, to work out the pattern of node numbers for the (i,j)
      % element. (See the solutions sheet description). With this worked
      % out, we can implement it as below.
      els(counter,:) = [(m+1)*(j-1)+i,...
                        (m+1)*(j-1)+i+1,...
                        (m+1)*j+i+1,...
                        (m+1)*j+i];
   end
end