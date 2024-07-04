%    """
%        ||> ------------------------------------------------------
%        ||>|                                                     |→→→→
%        ||>|                                                     |→→→        
%        ||>|                                                     |→→
%        ||>|                                                     |→
%        ||>|                                                    ←|
%        ||>|                                                   ←←|
%        ||>|                                                  ←←←|        
%        ||>|                                                 ←←←←|
%        ||> ------------------------------------------------------
%    """
%
% Author: YaLin Li, Department of Engineering Mechanics, 
% School of Civil Engineering, Southeast University
% tumu16lyl@163.com.

clear

plotres = 1; % Switch to control plotting of results

% Material parameters
E = 1000;
nu = 0.3;
conf = 1; % Plane strain

% Construct geometry
L = 0.5;
h = 0.1;
b = 1;
m = 200; n = 40;
[nodes,els] = quadMesh(m,n,L,h);
els = quad2triMesh(els); % Comment out to use quad elements, if available

% Find left face nodes
ln = find(nodes(:,1)==0);   % 找到x=0的坐标编号
ln1 = find(nodes(:,1)==0 & nodes(:,2)==0);


% Assemble BCs
bcInd1 = [2*ln-1; 2*ln1]'; % X- & Y-dofs of left face nodes  %找到ln中包含节点的1，2自由度;2*ln
bcVal1 = zeros(size(bcInd1)); % All of them fixed    %抓取该bcInd的矩阵，用Zeros定义成0矩阵

% 读取节点位移边界条件
nodaldisps = xlsread('/point_disp2.csv'); % 读取节点位移数据
dispNodeIndices = nodaldisps(:,1); % 获取节点号
dispX = nodaldisps(:,2); % 获取x方向位移
%dispY = nodaldisps(:,3); % 获取y方向位移

% Assemble BCs for displacement nodes
dispBcInd = [2*dispNodeIndices-1]'; % X- & Y-dofs of displacement nodes
dispBcVal = [dispX]'; % 位移值

% Combine boundary conditions
bcInd = [bcInd1 dispBcInd]; % 合并边界条件索引
bcVal = [bcVal1 dispBcVal]; % 合并边界条件值

% Point load
%pn = find(nodes(:,1)==L & nodes(:,2)==h/2); % Loaded node 找到x=L，y=h的节点编号
%plInd = 2*pn; % Just one node in this case (can be vect in general), y-dof affected 给该节点的y向自由度施加荷载
%plVal = -100; % Same here

% 线性分布载荷
%rightNodes = find(nodes(:,1) == L); % 找到x=L的所有节点
%numRightNodes = length(rightNodes);

% 从外部文件加载节点荷载

%nodalLoads = xlsread('/point_force.csv');      %节点力;

% 确保nodalLoads的尺寸正确（每行三个值：节点号，x方向荷载，y方向荷载）
%if size(nodalLoads, 2) ~= 3
%    error('nodalLoads文件的每行必须包含三个值：节点号，x方向荷载，y方向荷载');
%end

% 初始化载荷值
%plInd = zeros(2 * size(nodalLoads, 1)-1, 1);
%plVal = zeros(2 * size(nodalLoads, 1), 1);

% 填充载荷值
%for i = 1:size(nodalLoads, 1)
%    nodeId = nodalLoads(i, 1); % 获取节点号
%    plInd(2*i-1) = 2*nodeId-1; % x方向的载荷
%    plInd(2*i) = 2*nodeId;     % y方向的载荷
%    plVal(2*i-1) = nodalLoads(i, 2); % x方向的载荷值
%    plVal(2*i) = nodalLoads(i, 3);   % y方向的载荷值
%end

% Compile model
model.Mesh.Nodes = nodes;
model.Mesh.Elements = els;
model.Mesh.Thickness = b;
model.BoundaryConditions.Ind = bcInd;
model.BoundaryConditions.Val = bcVal;
%model.PointLoad.Ind = plInd;
%model.PointLoad.Val = plVal;
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
   
   % 创建一个新的figure
   figure;
   
   % 设置子图的位置和大小
   subplotPos = [...
       0.05 0.55 0.27 0.4;  % [left bottom width height]
       0.37 0.55 0.27 0.4;
       0.69 0.55 0.27 0.4;
       0.05 0.05 0.27 0.4;
       0.37 0.05 0.27 0.4;
       0.69 0.05 0.27 0.4];

   % X方向位移
   subplot('Position', subplotPos(1,:));
   patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',u(:,1), 'EdgeColor', 'none'); 
   axis equal;
   title('Disp x');
   colorbar; 
   colormap jet;
   caxis([-0.0228 0.0228]); % 设置色阶上下限 

   % Y方向位移
   subplot('Position', subplotPos(2,:));
   patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',u(:,2), 'EdgeColor', 'none'); 
   axis equal;
   title('Disp y');
   colorbar; 
   colormap jet;
   caxis([-0.1138 0]); % 设置色阶上下限

   % 位移大小
   %um = sqrt(u(:,1).*u(:,1) + u(:,2).*u(:,2));
   %subplot('Position', subplotPos(3,:));
   %patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',um, 'EdgeColor', 'none'); 
   %axis equal;
   %title('Disp mag');
   %colorbar; 
   %colormap jet;

   % 计算节点应力和应变
   uvec = reshape(u',size(nodes,1)*2,1);
   [En,Sn]=nodalStrainsStresses(nodes,els,uvec,E,nu,conf); % Nodal stresses and strains

   % X方向应力
   subplot('Position', subplotPos(3,:));
   patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',Sn(:,1), 'EdgeColor', 'none'); 
   axis equal;
   title('X Stress');
   colorbar; 
   colormap jet;
   caxis([-50 50]); % 设置色阶上下限

   % Y方向应力
   subplot('Position', subplotPos(4,:));
   patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',Sn(:,2), 'EdgeColor', 'none'); 
   axis equal;
   title('Y Stress');
   colorbar; 
   colormap jet;
   %caxis([-2 2]); % 设置色阶上下限

   % XY方向应力
   subplot('Position', subplotPos(5,:));
   patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',Sn(:,3), 'EdgeColor', 'none'); 
   axis equal;
   title('XY Stress');
   colorbar; 
   colormap jet;
   %caxis([-2 2]); % 设置色阶上下限  

   % Von Mises应力
   %figure;
   %Se = (0.5 * ( (Sn(:,1)-Sn(:,2)).^2+Sn(:,1).^2+Sn(:,2).^2+6*Sn(:,3).^2 ) ).^.5; % Von Mises stress
   %patch('Vertices',nodes+scalefactor*u,'Faces',els,'FaceColor','interp','FaceVertexCData',Se, 'EdgeColor', 'none'); 
   %axis equal;
   %title('Von Mises Stress');
   %colorbar; 
   %colormap jet;
   % 保存应力数据
   stress_data2 = table(nodes(:,1), nodes(:,2), Sn(:,1), Sn(:,2), Sn(:,3), ...
       'VariableNames', {'X', 'Y', 'Sigma_x', 'Sigma_y', 'Tau_xy'});
   writetable(stress_data2, 'stress_data2.csv');

   % 保存位移数据
   displacement_data2 = table(nodes(:,1), nodes(:,2), u(:,1), u(:,2), ...
       'VariableNames', {'X', 'Y', 'Disp_x', 'Disp_y'});
   writetable(displacement_data2, 'displacement_data2.csv');

end
