%SS-RFLE 使用模糊粗糙相似关系构造半监督信息权重（加入了信息熵构造的属性重要度sig）
function [ER,EC]=SSREFshu(H,f,k,d,q)
% INPUT
% data ------------- dataset(datasize*attributesize)
% f    ------------- number of cluster
% k    ------------- number of points in neighbors
% d    ------------- dimension
% q    ------------- rate of study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% clear all
% close all
% clc
%-----------------------------------------------------------
%----------------------------------------------------------
% f=2;       %类数
% q=0.5;     %降维学习速率
% k=5;      %邻域中点个数
% d=2;
% mu = [2 3 2 4 6];
% SIGMA = [1 0 0 0 0; 0 4 0 0 0;0 0 6 0 0;0 0 0 9 0;0 0 0 0 7];
% r1 = mvnrnd(mu,SIGMA,100);
% mu = [10 12 15 13 10];
% SIGMA = [2 0 0 0 0; 0 1 0 0 0;0 0 3 0 0;0 0 0 4 0;0 0 0 0 5];
% r2 = mvnrnd(mu,SIGMA,100);
% Z=[r1;r2];
% Zm=size(Z,1);
% l=[ones(5,1);zeros(Zm/2-5,1);2*ones(5,1);zeros(Zm/2-5,1)];
% H=[l Z];
%-------------------------------------------------------------
%-------------------------------------------------------------
% f=2;       %类数
% q=0.5;     %降维学习速率
% k=10;      %邻域中点个数
% d=3;
% Z=load('H:\叶东升硕士\eassy\Z.dat');
% [Zm,Zn]=size(Z);
% L=[ones(Zm/3,1);2*ones(Zm/3,1);zeros(Zm/3,1)];
% data0=xlsread('H:\叶东升硕士\eassy\data\seeds.xls');
% datar=data0(:,2:8);
% Z=datar;
% [Zm,Zn]=size(Z);
% l=[ones(10,1);zeros(140-20,1);2*ones(10,1);zeros(70,1)];
% H=[l Z];
%%
%X为带类标的数据集
%f为类数
%u为行向量，每类取得训练样本个数

[a,b]=size(H);
% [H0,H1,H]=fptesttrain(X,f,u);
%H0为训练集
%H1为测试集
%H为全体数据集，第一列为类标（为标记的样本的类标为0）
% H00=H0(:,2:b+1);

[P0,P1,HH,HH1]=fenlei(H,f);%其实是对训练集进行的fenlei
%%%%%%%%%%%%将数据X每个属性范围（0,1]
A=H(:,2:b);
A0=mean(A);
A1=std(A);
for i=1:a
    X1(i,2:b)=(A(i,:)-A0)./A1;
end
X1(:,1)=H(:,1);
% X1=H;
Z=X1(:,2:b);

%%%%%%%%%%%属性重要度的测量
l=P0';
labset1=find(l==1);
labset2=find(l==2);
ZCte=[Z(labset1,:);Z(labset2,:)];
labnum=[labset1;labset2];
HD=-length(labset1)/length(labnum)*log(length(labset1)/length(labnum))-length(labset2)/length(labnum)*log(length(labset2)/length(labnum));
HD=log(2);
for shui=1:b-1
    [aa,bb]=sort(ZCte(:,shui));
    for yangj=1:size(ZCte,1)-1
      c(1,1)=length(find(l(bb(1:yangj))==1));
      c(1,2)=yangj-c(1,1);
      c(2,1)=length(find(l(bb(1+yangj:end))==1));
      c(2,2)=size(ZCte,1)-yangj-c(2,1);
      caad=c./([yangj;size(ZCte,1)-yangj]*ones(1,2));
      for i=1:2
          for j=1:2
              if caad(i,j)==0;
                  logcaad(i,j)=0;
              else
                  logcaad(i,j)=log(caad(i,j));
              end
          end
      end
      GR(yangj,shui)=(HD+[yangj/size(ZCte,1),1-yangj/size(ZCte,1)]*sum(caad.*logcaad,2))/-(yangj/size(ZCte,1)*log(yangj/size(ZCte,1))+(size(ZCte,1)-yangj)/size(ZCte,1)*log((size(ZCte,1)-yangj)/size(ZCte,1)));
    end
end
G=max(GR);
SIG=G/sum(G)*(b-1);

%%%
% SIG=ones(1,b-1);
%%%%%%%%%%%利用属性重要度计算相似矩阵
R1=zeros(a,a);
R2=zeros(a,a);
for i=1:a
    for j=1:a
        for i0=2:b
            R1(i,j)=SIG(1,i0-1)*((X1(i,i0)-X1(j,i0))^2);
            R2(i,j)=R1(i,j)+R2(i,j);         
        end
     end
end
R2=sqrt(R2);
c1=1/max(max(R2));
R3=1-R2.*c1;%得到模糊相似矩阵
for i=1:a
   for j=1:a
       if X1(i,1)~=0&&X1(j,1)~=0
           if X1(i,1)==X1(j,1)
               R3(i,j)=1;
           else
               R3(i,j)=0;
           end
        end
    end
end
%%%%%%%%%%%%每一个样本属于标记的相同类的样本的隶属度应该是相同的
for i=1:f
   hh=HH1{i};
   R31=R3(:,hh);
   r31=(max(R31'))';%%%%取最大
   for j=1:a
        for j1=1:P1(i)
           R3(j,hh(1,j1))=r31(j,1);
        end
    end
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%利用属性重要度找k近邻
 X3=X1(:,2:b);           %利用归一化的数据 
 D=zeros(a,a);
    for i=1:a
        D0=(repmat(X3(i,:),[a,1])-X3)';
        D(i,:)=SIG*(D0.^2);
    end
    D=sqrt(D);
    [D1,D2]=sort(D);
    N=D2(2:(k+1),:);%以上构造k近邻
    N=N';
    %%%%%%%%%%%%%找近邻下近似模糊粗糙集UN
    UN=zeros(a,a);
    for i=1:a
        for j=1:a
             UN(j,i)=min(R3(N(j,:),i));%/length(N(j,:));%%%%%%%%%%%%%%%%%
        end
    end
    UN;
    %%%%%%%%%%%%%近邻赋权值
    sigma=(sum(sum(D.^2))) /(a^2);
    WN1=exp(-(D.^2)/sigma);
    WN2=zeros(a,a);
    %-------------------------------------论文
%     for i=1:a
%         for j=1:k
%             if H(i,1)==0
%                 WN2(i,N(i,j))=(WN1(i,N(i,j))*UN(i,N(i,j)));
%             else
%                 if H(N(i,j),1)==1
%                      WN2(i,N(i,j))=(WN1(i,N(i,j))*UN(i,N(i,j)));
%                 elseif H(i,1)==H(N(i,j),1)
%                      WN2(i,N(i,j))=1;
%                 else
%                      WN2(i,N(i,j))=0;
%                 end
%             end
%          end
%     end
%     WN2;
%     WN0=eye(a,a);
%     WN=max(WN2,WN2');
%     WN=WN+WN0;
    %--------------------------------------
    %--------------------------------------eassy
    for i=1:f
       hh=HH1{i};
       Nf{i}=[];
       for j=1:length(hh)
       Nf{i}=[Nf{i};N(hh,:)];
       end
    end
     for i=1:f
         INf{i}=[];
         for j=1:f
             if j==i
                  NN{j}=[];
             else
                  NN{j}=Nf{j};
             end
        INf{i}=[INf{i};NN{j}];
         end
     end
    for i=1:a
        for j=1:a
            if H(i,1)~=0&&H(j,1)~=0&&H(i,1)==H(j,1)
                WN2(i,j)=1;
            else if H(i,1)~=0&&H(j,1)~=0&&H(i,1)~=H(j,1)
                    WN2(i,j)=0;
                else if ismember(i,N(j,:))||ismember(j,N(i,:))
                        WN2(i,j)=max(UN(i,j),UN(j,i))*WN1(i,j);
                        else if (H(i,1)~=0&&ismember(j,Nf{P0(i)})&&~ismember(j,INf{P0(i)}))||(H(j,1)~=0&&ismember(i,Nf{P0(j)})&&~ismember(i,INf{P0(j)}))
                                WN2(i,j)=max(UN(i,j),UN(j,i))*WN1(i,j);
                            else
                                WN2(i,j)=0;
                            end
                    end
                end
            end
        end
    end
  WN=WN2;
    %--------------------------------------
    %%%%%%%%%%%%所有样本到每类的隶属度矩阵
    WC=zeros(f,a);
    for i1=1:f
        XC=[];
        XC=HH1{i1};
        WC(i1,:)=R3(:,XC(1,1));%因为相同标记的在R3中是相同的，取每类中一个即可
    end
    WC;

    %%%%%%%%%%%%%%%%
    W=zeros(a+f,a+f);
    W(1:f,1:f)=eye(f);
    W(1:f,f+1:f+a)=q*WC;
    W(f+1:f+a,1:f)=q*WC';
    W(f+1:f+a,f+1:f+a)=2*(1-q)*WN;

    D=diag(sum(W,2));%W为权矩阵
    L=D-W;
    [ER0,EV0]=eig(L,D);%%LY'=rDY'
    EV1=diag(EV0);


    [EV11,EV12]=sort(EV1,'ascend');
    
    t=1;
    while t<a   
        if EV11(t)>1e-5 %%去非零值
             ER1=zeros(a+f,d);%%d个特征值对应的d个特征向量
             for i=t:d+t-1
                  ER1(:,i)=ER0(:,EV12(i));
             end
             ER2=ER1(:,t:d+t-1);
             break
        else 
             t=t+1;
        end
    end
    ER2;
    ER=ER2(f+1:a+f,:); %数据在低维的表示
    EC=ER2(1:f,:);  %数据集在低维的类中心