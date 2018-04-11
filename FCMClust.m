function [center, U, obj_fcm] = FCMClust(data,cluster_n,options)
% FCMClust.m   采用模糊C均值对数据集data聚为cluster_n类 
% 用法：
%   1.  [center,U,obj_fcm] = FCMClust(Data,N_cluster,options);
%   2.  [center,U,obj_fcm] = FCMClust(Data,N_cluster);
% 输入：
%   data        ---- nxm矩阵,表示n个样本,每个样本具有m的维特征值
%   N_cluster   ---- 标量,表示聚合中心数目,即类别数
%   options     ---- 4x1矩阵，其中
%       options{1}   :  初始聚类中心                          （缺省：随机隶属度矩阵）
%       options{2}(1):  隶属度矩阵U的指数，>1                  (缺省值: 2.0)
%       options{2}(2):  最大迭代次数                           (缺省值: 100)
%       options{2}(3):  隶属度最小变化量,迭代终止条件           (缺省值: 1e-5)
%       options{2}(4):  每次迭代是否输出信息标志                (缺省值: 1)
% 输出：
%   center      ---- 聚类中心
%   U           ---- 隶属度矩阵
%   obj_fcm     ---- 目标函数值
%   Example:
%       data = rand(100,2);
%       [center,U,obj_fcm] = FCMClust(data,2);
%       plot(data(:,1), data(:,2),'o');
%       hold on;
%       maxU = max(U);
%       index1 = find(U(1,:) == maxU);
%       index2 = find(U(2,:) == maxU);
%       line(data(index1,1),data(index1,2),'marker','*','color','g');
%       line(data(index2,1),data(index2,2),'marker','*','color','r');
%       plot([center([1 2],1)],[center([1 2],2)],'*','color','k')
%       hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_n = size(data, 1); % 求出data的第一维(rows)数,即样本个数
in_n = size(data, 2);   % 求出data的第二维(columns)数，即特征值长度

% 默认操作参数
% default_options{1}=options{1};
default_options = [2;	              % 隶属度矩阵U的指数
                      100;                % 最大迭代次数 
                      1e-5;               % 隶属度最小变化量,迭代终止条件
                      1];                 % 每次迭代是否输出信息标志
if  length(options{2})~=4;
	options{2} = default_options;
end     %分析有options做参数时候的情况
%将options 中的分量分别赋值给四个变量;
expo = options{2}(1);          % 隶属度矩阵U的指数
max_iter = options{2}(2);		% 最大迭代次数 
min_impro = options{2}(3);		% 隶属度最小变化量,迭代终止条件
display = options{2}(4);		% 每次迭代是否输出信息标志 
obj_fcm = zeros(max_iter, 1);	% 初始化输出参数obj_fcm
% if size(U,2)~=size(data,1)
if size(options{1},1)==0
    U = initfcm(cluster_n, data_n);     % 初始化模糊分配矩阵,使U满足列上相加为1,
else dist = distfcm(options{1}, data);       % 计算距离矩阵
    tmp = dist.^(-2/(expo-1));     
    U= tmp./(ones(cluster_n, 1)*sum(tmp));
end

% Main loop  主要循环
for i = 1:max_iter
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
	[U, center, obj_fcm(i)] = stepfcm(data, U, cluster_n, expo);
	if display
		fprintf('FCM:Iteration count = %d, obj. fcm = %f\n', i, obj_fcm(i));
	end
	% 终止条件判别
	if i > 1
		if abs(obj_fcm(i) - obj_fcm(i-1)) < min_impro
            break;
        end
	end
end

iter_n = i;	% 实际迭代次数 
obj_fcm(iter_n+1:max_iter) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 子函数1
% 
function U = initfcm(cluster_n, data_n)
% 初始化fcm的隶属度函数矩阵
% 输入:
%   cluster_n   ---- 聚类中心个数
%   data_n      ---- 样本点数
% 输出：
%   U           ---- 初始化的隶属度矩阵
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 子函数2
% 
function [U_new, center, obj_fcm] = stepfcm(data, U, cluster_n, expo)
% 模糊C均值聚类时迭代的一步
% 输入：
%   data        ---- nxm矩阵,表示n个样本,每个样本具有m的维特征值
%   U           ---- 隶属度矩阵
%   cluster_n   ---- 标量,表示聚合中心数目,即类别数
%   expo        ---- 隶属度矩阵U的指数                      
% 输出：
%   U_new       ---- 迭代计算出的新的隶属度矩阵
%   center      ---- 迭代计算出的新的聚类中心
%   obj_fcm     ---- 目标函数值
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))'); % 新聚类中心(5.4)式
dist = distfcm(center, data);       % 计算距离矩阵
obj_fcm = sum(sum((dist.^2).*mf));  % 计算目标函数值 (5.1)式
tmp = dist.^(-2/(expo-1));     
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % 计算新的隶属度矩阵 (5.3)式
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 子函数3
% 
function out = distfcm(center, data)
% 计算样本点距离聚类中心的距离
% 输入：
%   center     ---- 聚类中心
%   data       ---- 样本点
% 输出：
%   out        ---- 距离
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1), % 对每一个聚类中心
    % 每一次循环求得所有样本点到一个聚类中心的距离
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
%    out(k, :) = max(abs(data-ones(size(data,1),1)*center(k,:))');
end