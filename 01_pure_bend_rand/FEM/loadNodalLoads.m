function nodalLoads = loadNodalLoads(filename)
% loadNodalLoads 读取CSV文件并返回节点荷载向量
%   filename: 输入文件名
%
%   nodalLoads: 返回的节点荷载向量，大小为 2*nN-by-1

    % 读取CSV文件
    data = readmatrix(filename);
    
    % 获取节点数
    nodeIndices = data(:,1);
    forcesX = data(:,2);
    forcesY = data(:,3);
    
    % 获取所有节点的数量
    nN = max(nodeIndices);
    
    % 初始化荷载向量
    nodalLoads = zeros(2*nN, 1);
    
    % 遍历每个节点并分配荷载
    for i = 1:length(nodeIndices)
        nodeIndex = nodeIndices(i);
        nodalLoads(2*nodeIndex-1) = forcesX(i); % x方向荷载
        nodalLoads(2*nodeIndex) = forcesY(i);   % y方向荷载
    end
end
