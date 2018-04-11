function [center, U, dis] = SSPCMClustC(datar,l,cluster_n,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% clc
consi=0;
% % mu = [2 3];
% % SIGMA = [1 0;0 10];
% % r1 = mvnrnd(mu,SIGMA,100);
% % mu = [10 12];
% % SIGMA = [2 0;0 12];
% % r2 = mvnrnd(mu,SIGMA,100);
% % scatter(r1(:,1),r1(:,2),'r'),hold on
% % scatter(r2(:,1),r2(:,2),'b')
% % data=[r1;r2] ;
%%%%%%%%%%%%%%%%%%%%%%
%SHIYAN
% Rd=load('H:\Rd.dat');
% C=load('H:\C.dat');
% l=load('H:\l.dat');
% datar=Rd;
% cluster_n=2;
% options={[C],[]};
%%%%%%%%%%%%%%%%%%%%%%
% cluster_n=2;
% data0=xlsread('seeds.xls');
% datar=data0(:,2:8);
% Lr=data0(:,1);
% dat1=find(Lr==1);
% dat2=find(Lr==2);
% dat3=find(Lr==3);
% % data=[data(:,4) data(:,1:3)];
% % scatter(datar(:,1),datar(:,2),'r')
if consi==1
[data,PS]=mapstd(datar');
data=data';
else
    data=datar;
end
[data_n,in_n]=size(data);
% l=[ones(50,1);zeros(data_n-100,1);2*ones(50,1)];
% % scatter(data(:,1),data(:,2),'b')
% % figure()
% for i=1:6
%     for j=i+1:7
% scatter(data(dat1,i),data(dat1,j),'r'),hold on
% scatter(data(dat2,i),data(dat2,j),'g'),hold on
% scatter(data(dat3,i),data(dat3,j),'b')
%  title(['i= ',num2str(i), ' , j= ',num2str(j)])
% pause(0.01)
% hold off
%     end
% end
% xlswrite('E:\eassy\data.xls',data)
% data0=xlsread('E:\eassy\data.xls');
% r3= [2*rand(100,1);3*rand(100,1)];
% data=[data0,r3];
% l=[zeros(70,1);ones(10,1);zeros(60,1);2*ones(10,1);zeros(60,1)];
% l=[zeros(10,1);zeros(100-20,1);1*ones(10,1);2*ones(10,1);zeros(40,1)];
% l=[ones(10,1);zeros(140-20,1);2*ones(10,1);zeros(70,1)];
% options={rand(cluster_n,in_n),[]};
K=1;


%---------------------------------------------
default_options    = [2;	                  % 隶属度矩阵U的指数
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
center=options{1};
U=rand(cluster_n, data_n);
W=1/in_n*ones(in_n,1);
%---------------------------------------------------
%---------------------------------------------------
bata1=sum(sum((data-repmat(sum(data)/data_n,[data_n,1])).^2))/data_n;
bata2=bata1;
% bata2=sum(sum(data.^2));
a=zeros(cluster_n,data_n);
 for ii=1:data_n
  if l(ii)~=0
   a(:,ii)=ones(cluster_n,1);
   a(l(ii),ii)=-1;
  end
 end
i=2;
while i<=max_iter
   
    %在第k步循环中改变聚类中心ceneter,和分配函数U的隶属度值;
mf = U.^expo;       % 隶属度矩阵进行指数运算结果
% for k = 1:cluster_n, % 对每一个聚类中心
%     dist(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
% end    % 计算距离矩阵
% ata=sum(mf.*dist.^2,2)./sum(mf,2);      %计算系数ata
% obj_fcm(i) = sum(sum((dist.^2).*mf))+sum(ata.*sum((ones(cluster_n,data_n)-U).^expo,2));  % 计算目标函数值  
% U =ones(cluster_n,data_n)./(ones(cluster_n,data_n)+(dist.^2./repmat(ata,[1,data_n])).^(1/(expo-1)));   % 计算新的隶属度矩阵
Vr=zeros(in_n,cluster_n);
for k=1:cluster_n
dif{k}=(data-ones(data_n,1)*center(k,:)).^2;
bbb(:,k)=sum( mf(k,:)'*(1-W').*dif{k}/bata2)';
% pause
Vr(:,k)=exp(-sum(mf(k,:)'*(1-W').*dif{k}/bata2))';
end
V=Vr./repmat(sum(Vr),[in_n,1]);
Wr=zeros(in_n,cluster_n);
for k=1:cluster_n
    Wr(:,k)=(-sum(mf(k,:)'*(1-V(:,k)').*dif{k}))';
end
W=exp(sum(Wr,2)/bata2/cluster_n)/sum(exp(sum(Wr,2)/bata2/cluster_n));
dis=zeros(cluster_n,data_n);
for k=1:cluster_n
    dis(k,:)=((in_n/max((in_n-2),1))^2*dif{k}*((1-V(:,k)).*(1-W)))';
end

Ur=U;
U=exp(-expo*cluster_n^1/2/bata1*(1+K*a).*dis);
center = (U.^expo*data)./repmat(sum(U.^expo,2),[1,in_n]); % 新聚类中心(5.4)式


% a1=sum(sum(dis.*mf)),
% a2=bata1/expo^2/cluster_n^(1/2)*sum(sum(mf.*log(mf)-mf)),
% a3=sum(sum(K*dis.*a.*mf)),
% a4=-bata2*sum(V.*log(V),2),
% a5=-bata2*(W.*log(W))
obj_fcm(i)=sum(sum(dis.*mf))+bata1/expo^2/cluster_n^(1/2)*sum(sum(mf.*log(mf)-mf))+sum(sum(K*dis.*a.*mf))-bata2*sum(sum(V.*log(V)))-cluster_n*bata2*(W'*log(W));
caa=max(max(abs(U-Ur)));
	if display
		fprintf('PCM:Iteration count = %d, obj. pcm = %f\n', i, obj_fcm(i));
	end
	% 终止条件判别

% 		if abs(obj_fcm(i) - obj_fcm(i-1)) < min_impro
        if caa < min_impro
            break;      
        end
    i=i+1;
end
% if consi==1
% datac=[data;center];
% datacr=(mapstd('reverse',datac',PS))';
% center=datacr(end-cluster_n+1:end,:);
% end
iter_n = i;	% 实际迭代次数
obj_fcm(iter_n+1:max_iter) = [];
% Lc=zeros(length(l),1);
% Cr1=find(U(1,:)>=U(2,:));
% for ii=1:length(Cr1)
%     if U(1,Cr1(ii))<0.1
%         Lc(Cr1(ii))=3;
%     else if U(2,Cr1(ii))>0.8
%           [zz1,aa1]=min([dis(1,Cr1(ii)),dis(2,Cr1(ii))]);
%           Lc(Cr1(ii))=aa1;
%         end
%     end
% end
% Cr2=find(U(1,:)<U(2,:));
% for jj=1:length(Cr2)
%     if U(2,Cr2(jj))<0.1
%         Lc(Cr2(jj))=3;
%     else if U(1,Cr2(jj))>0.8
%           [zz2,aa2]=min([dis(1,Cr2(jj)),dis(2,Cr2(jj))]);
%           Lc(Cr2(jj))=aa2;
%         end
%     end
% end 
%-----------------------------------------------------------


%---------------------------------------------------------

% class(1)=length(find(Lc(dat1)==1));
% class(2)=length(find(Lc(dat2)==2));
% class(3)=length(find(Lc(dat3)==3));
%---------------------
%hutu
% d1=find(Lc==1);
% d2=find(Lc==2);
% d3=find(Lc==3);
% for i=1:in_n-1
%     for j=i+1:in_n
%         figure()
% scatter(data(d1,i),data(d1,j),'r'),hold on
% scatter(data(d2,i),data(d2,j),'g'),hold on
% scatter(data(d3,i),data(d3,j),'b')
%  title(['i= ',num2str(i), ' , j= ',num2str(j)])
% pause(0.01)
% scatter(center(1,i),center(1,j),'k+')
% scatter(center(2,i),center(2,j),'k+')
% hold off
%     end
% end
%--------------------