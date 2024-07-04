function [KK,FF] = imposeBCs(K,F,BC_ind,BC_val)
% Adjust the global stiffness matrix K and force vector F to reflect a set
% of boundary conditions for a particular problem. For an n-DOF problem, K
% is n-by-n, F is n-by-1. KK and FF are the modified versions of K and F.
%
% BC_ind is a list of degrees of freedom to be constrained, and BC_val is a
% list of corresponding values for those degrees of freedom. For example if
% we have a 2D stress analysis problem and we want to fix the
% x-displacements of nodes 5 through to 8, we would provide the following:
% BC_ind = [9 11 13 15]
% BC_val = [0 0 0 0]
%
% If DOF i is to be constrained, the stiffness matrix K is modified so that
% row i and column i are all zeroes, except for entry (i,i), which retains
% its value. The force vector F is then also modified to reflect these
% changes.
% 
% A four step process is followed:
% 1. For each BC, the corresponding rows are zeroed, except for the
% diagonal entries.
% 2. The corresponding force values are updated accordingly.
% 3. All other force values are updated accordingly.
% 4. Corresponding columns in K are zeroed, again, except for the diagonal
% entries.
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, April 2013.

KK = K;
FF = F;
% Adjust rows
for i = 1:length(BC_ind)
   KK(BC_ind(i),:) = 0;
   KK(BC_ind(i),BC_ind(i)) = K(BC_ind(i),BC_ind(i));
end
% Adjust corresponding force components
for i = 1:length(BC_ind)
   FF(BC_ind(i)) = KK(BC_ind(i),BC_ind(i)) * BC_val(i);
end
% Adjust remaining force components
for i = 1:length(F)
   if ~any(BC_ind == i) % only want to adjust components NOT corresponding to constrained DOFs
      FF(i) = F(i) - BC_val*K(i,BC_ind)';
   end
end
% Adjust cols
for i = 1:length(BC_ind)
   KK(:,BC_ind(i)) = 0;
   KK(BC_ind(i),BC_ind(i)) = K(BC_ind(i),BC_ind(i));
end