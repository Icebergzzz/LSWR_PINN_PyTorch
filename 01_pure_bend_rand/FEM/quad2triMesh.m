function tE = quad2triMesh(qE)
% Convert a mesh of quadrilateral elements defined by element array qE into
% an equivalent mesh of triangular elements tE.
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, Feb 2014.

n = size(qE,1);
tE = zeros(n,3);
for i = 1:n
   % Each quad is split into two tris. If the quad nodes are numbered
   % [i,j,k,l], the two tris will be [i,j,k] & [k,l,i]
   tE(2*i-1,:) = [qE(i,1) qE(i,2) qE(i,3)]; % First tri
   tE(2*i,:) = [qE(i,3) qE(i,4) qE(i,1)]; % Second tri
end


%function tE = quad2triMesh(qE)
% Convert a mesh of quadrilateral elements defined by element array qE into
% an equivalent mesh of quadrilateral elements tE.
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, Feb 2014.

%n = size(qE,1);
%tE = zeros(n,4);
%for i = 1:n
   % Each quad is split into two quads. If the quad nodes are numbered
   % [i,j,k,l], the two quads will be [i,j,k,m] & [k,l,m,i]
   %tE(2*i-1,:) = [qE(i,1) qE(i,2) qE(i,3) qE(i,4)]; % First quad
   %tE(2*i,:) = [qE(i,3) qE(i,4) qE(i,1)]; % Second quad
%end