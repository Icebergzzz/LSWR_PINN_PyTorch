% An example FE problem solved using fem2d. Key points demonstrated:
% 1. Mesh generation
% 2. Boundary condition construction
% 3. Load construction
% 4. Constructing the 'model' struct for input to fem2d
%
% Problem description:
% Cantilever beam fixed on left hand end, with a point load applied at the
% free end.
% Length L = 1 m, rectangular cross-section with dimensions b = 100 mm, h =
% 60 mm.
% Material properties: E = 1000 Pa, nu = 0.3.
% Load: P = 1000 N.
% Plane stress condition assumed.
%
%                                           P
%                                           |
%  / _______________________________________V
%  /|                                       |
%  /|                                       |
%  /|_______________________________________|                                      
%  /
%
% Author: Zeike Taylor, Dept of Mechanical Engineering, University of
% Sheffield, z.a.taylor@sheffield.ac.uk, January 2017.

clear

plotres = 1; % Switch to control plotting of results

% Material parameters
E = 1000;
nu = 0.3;
conf = 1; % Plane strain

% Construct geometry
L = 0.5;
h = 0.1;
%D = h;
b = 1   ;
m = 100; n = 20;
[nodes,els] = quadMesh(m,n,L,h);
els = quad2triMesh(els); % Comment out to use quad elements, if available

% Find left face nodes
ln = find(nodes(:,1)==0);   % 找到x=0的坐标编号

% Assemble BCs
bcInd = [2*ln-1;2*ln]'; % X- & Y-dofs of left face nodes  %找到ln中包含节点的1，2自由度
bcVal = zeros(size(bcInd)); % All of them fixed    %抓取该bcInd的矩阵，用Zeros定义成0矩阵

% 线性分布载荷：顶部 -50 N，底部 +50 N
rightNodes = find(nodes(:,1) == L); % 找到x=L的所有节点
numRightNodes = length(rightNodes);
%loadTop = -50 * 0.001; % 顶部载荷
%loadBottom = 50 * 0.001; % 底部载荷
%stressTop = -50; % 顶部应力
%stressBottom = 50; % 底部应力

% 从外部文件加载节点荷载

nodalLoads = xlsread('/point_force.csv');      %节点力;

% 确保nodalLoads的尺寸正确（每行三个值：节点号，x方向荷载，y方向荷载）
if size(nodalLoads, 2) ~= 3
    error('nodalLoads文件的每行必须包含三个值：节点号，x方向荷载，y方向荷载');
end

% 初始化载荷值
plInd = zeros(2 * size(nodalLoads, 1)-1, 1);
plVal = zeros(2 * size(nodalLoads, 1), 1);

% 填充载荷值
for i = 1:size(nodalLoads, 1)
    nodeId = nodalLoads(i, 1); % 获取节点号
    plInd(2*i-1) = 2*nodeId-1; % x方向的载荷
    plInd(2*i) = 2*nodeId;     % y方向的载荷
    plVal(2*i-1) = nodalLoads(i, 2); % x方向的载荷值
    plVal(2*i) = nodalLoads(i, 3);   % y方向的载荷值
end



%plInd = 2*rightNodes - 1; % 作用于每个节点的x自由度
%plVal = zeros(numRightNodes, 1); % 初始化每个节点的载荷值

% 遍历每个节点并将荷载添加到全局荷载向量中
%for i = 1:numRightNodes
%    plVal(2*i-1) = nodalLoads(2*i-1);
%    plVal(2*i) = nodalLoads(2*i);
%end


% Construct loaded surface (top half of left face)
%ln = find(nodes(:,1)==L);
%for i = 1:length(ln)-1
%   S(i,:) = [ln(i),ln(i+1)]; % Connectivity of line els along this face
%end

% 计算每个节点的载荷值
%for i = 1:numRightNodes
%    y = nodes(rightNodes(i), 2); % 获取节点的y坐标
%    plVal(i) =  (loadBottom - loadTop) * (y / h); % 线性插值载荷值
%end

% 计算每个节点的力值
%for i = 1:numRightNodes
%    y = nodes(rightNodes(i), 2); % 获取节点的y坐标
%    stress = (stressBottom - stressTop) * (y / h); % 线性插值应力值

    % 确定节点间的距离
%    if i == 1
%        dy = nodes(rightNodes(i+1), 2) - y; % 第一个节点到下一个节点的距离
%    elseif i == numRightNodes
%        dy = y - nodes(rightNodes(i-1), 2); % 最后一个节点到前一个节点的距离
%    else
%        dy = (nodes(rightNodes(i+1), 2) - nodes(rightNodes(i-1), 2)) / 2; % 中间节点的距离
%    end

%    plVal(i) = stress * b * dy; % 将应力转换为力，施加到节点
%end

% 计算每个节点的力值
%for i = 1:numRightNodes
%    y = nodes(rightNodes(i), 2); % 获取节点的y坐标
%    stress =  (stressBottom - stressTop) * (y / h); % 线性插值应力值
%    plVal(i) = stress * b * (h / numRightNodes); % 将应力转换为力，施加到节点
%end

% Pressure magnitude
%P_top = [-50,0];
%P_bottom = [50,0];

% Compile model
model.Mesh.Nodes = nodes;
model.Mesh.Elements = els;
model.Mesh.Thickness = b;
model.BoundaryConditions.Ind = bcInd;
model.BoundaryConditions.Val = bcVal;
model.PointLoad.Ind = plInd;
model.PointLoad.Val = plVal;


model.ConstitutiveModel.Young = E;
model.ConstitutiveModel.Poisson = nu;
model.ConstitutiveModel.Configuration = conf;

% Solve the model
u = fem2d(model);

% Postprocess - showing a mix of results in this case
% (These are just examples! For your project, it's up to you to decide what
%  is required.)
if plotres == 1
   scalefactor = 0; % To accentuate the deformations
   %figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','c'); axis equal
   %title('Deformed mesh')
   figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',u(:,1), 'EdgeColor', 'none'); axis equal
   title('Disp x'), colorbar, colormap jet
   caxis([-0.0228 0.0228]); % 设置色阶上下限 

   figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',u(:,2), 'EdgeColor', 'none'); axis equal
   title('Disp y'), colorbar, colormap jet
   caxis([-0.1138 0]); % 设置色阶上下限

   um = sqrt(u(:,1).*u(:,1) + u(:,2).*u(:,2));
   %figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',um); axis equal
   %title('Disp mag'), colorbar, colormap jet
   uvec = reshape(u',size(nodes,1)*2,1);
   [En,Sn]=nodalStrainsStresses(nodes,els,uvec,E,nu,conf); % Nodal stresses and strains

   %figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',En(:,1)); axis equal
   %title('X Strain'), colorbar, colormap jet

   figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',Sn(:,1), 'EdgeColor', 'none'); axis equal
   title('X Stress'), colorbar, colormap jet
   caxis([-50 50]); % 设置色阶上下限

   Se = (0.5 * ( (Sn(:,1)-Sn(:,2)).^2+Sn(:,1).^2+Sn(:,2).^2+6*Sn(:,3).^2 ) ).^.5; % Von Mises stress
   %figure, patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',Se); axis equal
   %title('Von Mises Stress'), colorbar, colormap jet

end