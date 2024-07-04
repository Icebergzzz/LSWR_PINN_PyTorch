function [nodes,els] = quarterAnnulusQ4Mesh(m,n,R1,R2)
% Construct a mesh of quadrilateral elements for a quarter annulus.
%
% Inputs:
% - m:  number of elements in circumferential direction
% - n:  number of elements in radial direction
% - R1: inner radius
% - R2: outer radius
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, October 2014.

% Create nodes
nodes = zeros((m+1)*(n+1),2);
dr = (R2-R1)/n; %  Radial node spacing
dt = pi/1/m; % Circumferential node angular spacing
cnt = 0; % Node counter
for i = 1:n+1 % Loop over rows (radial position)
   r = R1 + (i-1)*dr; % Radial position of current row
   for j = 1:m+1 % Loop over columns (angular position)
      cnt = cnt+1;
      t = (j-1)*dt; % Angular position of current column
      nodes(cnt,:) = [r*cos(t),r*sin(t)];
   end
end

% Create elements
els = zeros(m*n,4);
cnt = 0; % Now, counter keeps track of the element ids
for j = 1:n
   for i = 1:m
      cnt = cnt + 1;
      els(cnt,:) = [(m+1)*(j-1)+i,(m+1)*(j-1)+i+1,(m+1)*j+i+1,(m+1)*j+i];
   end
end
