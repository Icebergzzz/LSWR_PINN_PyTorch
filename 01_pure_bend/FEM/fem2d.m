function u = fem2d(model)
% Solve the 2D linear elastic finite element problem defined in model.
%
% model is a struct containing:
% - Mesh -> struct containing arrays Nodes, Elements and Thickness
% - PressureLoad (optional) -> struct containing arrays P and S
% - PointLoad (optional) -> struct containing arrays Ind and Val
% - BoundaryConditions -> struct containing arrays Ind and Val
% - ConstitutiveModel -> struct containing Young's modulus, Poisson's
%   ratio and configuration (plane strain or plane stress)
% 
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2017.

% ==================== Process inputs ====================

% Mesh
nodes = model.Mesh.Nodes;
els = model.Mesh.Elements;
t = model.Mesh.Thickness;
nnds = size(nodes,1);

% Point loads
if isfield(model,'PointLoad')
   plInd = model.PointLoad.Ind;
   plVal = model.PointLoad.Val;
end

% Pressure loads
if isfield(model,'PressureLoad')
   P = model.PressureLoad.P;
   S = model.PressureLoad.S;
end

% linearVaryingPressureLoad
if isfield(model,'linearVaryingPressureLoad')
   P_top = model.linearVaryingPressureLoad.P_top;
   P_bottom = model.linearVaryingPressureLoad.P_bottom;
   S = model.linearVaryingPressureLoad.S;
end

% BCs
bcInd = model.BoundaryConditions.Ind;
bcVal = model.BoundaryConditions.Val;

% Material props
E = model.ConstitutiveModel.Young;
nu = model.ConstitutiveModel.Poisson;
conf = model.ConstitutiveModel.Configuration;

% ==================== Construct system matrices ====================

% Stiffness matrix
K = stiffnessGlobal(nodes,els,t,E,nu,conf);

% External loads (only pressure loads supported in this code)
F = zeros(nnds*2,1);
if isfield(model,'PointLoad')
   F = F + pointLoad(plInd,plVal,nnds);
end
if isfield(model,'PressureLoad')
   F = F + uniformPressureLoad(P,S,nodes,t);
end
if isfield(model,'linearVaryingPressureLoad')
   F = F + linearVaryingPressureLoad(P_top, P_bottom, S, nodes, D, t);
end

% Apply BCs
[KK,FF] = imposeBCs(K,F,bcInd,bcVal);

% ==================== Solve ====================
UU = KK\FF;

u = reshape(UU,2,nnds)';
