function [center, U, obj_pcm] = PCMClustC(data,cluster_n,options)
% INPUT
% data ------------- dataset(datasize*attributesize)
% cluster_n----------- number of cluster
% options     ---- cell（1*2），
%       options{1}   :  initial center of cluster                    
%       options{2}(1):  exponent for the matrix U             (default: 2.0)
%       options{2}(2):  maximum number of iterations          (default: 100)
%       options{2}(3):  minimum amount of improvement         (default: 1e-5)
%       options{2}(4):  info display during iteration         (default: 1)
% OUTPUT
% center-------------  center of cluster
% U------------------  probability matrix
% obj_pcm------------  object function value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% clc
% mu = [2 3];
% SIGMA = [1 0;0 7];
% r1 = mvnrnd(mu,SIGMA,100);
% mu = [10 15];
% SIGMA = [2 0;0 5];
% r2 = mvnrnd(mu,SIGMA,100);
% scatter(r1(:,1),r1(:,2),'r'),hold on
% scatter(r2(:,1),r2(:,2),'b')
% data=[r1;r2];
% options={[2 3;10 15],[]};
% %----------------------------
% cluster_n=2;
% options=[2 100 1e-5 1];
% expo = options(1);          % 可能度矩阵U的指数
% max_iter = options(2);		% 最大迭代次数 
% min_impro = options(3);		% 可能度最小变化量,迭代终止条件
% display = options(4);		% 每次迭代是否输出信息标志 
% obj_fcm = zeros(max_iter+1, 1);	% 初始化输出参数obj_fcm
% [data_n,in_n]=size(data);
%  U = 0.5*ones(cluster_n, data_n);
%  center=[2 3;10 15];
% %----------------------------
%---------------------------------------------
%---------------------------------------------
[data_n,in_n]=size(data);
default_options = [2;	                  % 隶属度矩阵U的指数
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
obj_pcm = zeros(max_iter, 1);	% 初始化输出参数obj_fcm
center=options{1};
U=0.5*ones(cluster_n, data_n);
%---------------------------------------------------
%---------------------------------------------------
i=2;
while i<=max_iter
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
for k = 1:cluster_n, % 对每一个聚类中心
    dist(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end    % 计算距离矩阵
ata=sum(mf.*dist.^2,2)./sum(mf,2);      %计算系数ata
obj_pcm(i) = sum(sum((dist.^2).*mf))+sum(ata.*sum((ones(cluster_n,data_n)-U).^expo,2));  % 计算目标函数值  
U =ones(cluster_n,data_n)./(ones(cluster_n,data_n)+(dist.^2./repmat(ata,[1,data_n])).^(1/(expo-1)));   % 计算新的可能度矩阵
center = (U.^expo*data)./repmat(sum(U.^expo,2),[1,in_n]); % 新聚类中心(5.4)式
	if display
		fprintf('PCM:Iteration count = %d, obj. pcm = %f\n', i, obj_pcm(i));
	end
	% 终止条件判别

		if abs(obj_pcm(i) - obj_pcm(i-1)) < min_impro
            break;      
        end
    i=i+1;
%    pause(0.001)

end

iter_n = i;	% 实际迭代次数
obj_pcm(iter_n+1:max_iter) = [];
