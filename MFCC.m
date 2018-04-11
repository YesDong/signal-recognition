% close all
% clear all
% clc
% wav= audioread('work50-1-all 00_00_00-00_00_30.wav',[1 307712]);
function Z=cood1(wav,dim)
% eg. 200个文件 i_files= 1:200
% 307712是提取的样本数（限制读取的长度，约28s）
% 要提取的MFCC系数个数
wav=wav(:,dim);

num_ceps_coeffs = 14;
c.fs = 11025;  %采样频率
% 设置每帧大小（包含样本数）
c.seg_size = 1024; 
c.hop_size = 512; %% c.seg_size-交叠部分=c.hop_size 
% 帧数计算
num_segments = floor((length(wav)-c.seg_size)/c.hop_size)+1;           
% 初始化功率谱矩阵
P = zeros(c.seg_size/2+1,num_segments); 
% 设置窗函数
c.w = 0.5*(1-cos(2*pi*(0:c.seg_size-1)/(c.seg_size-1)))';%汉宁窗函数
%size(c.w)
% 逐帧做FFT            
for i_p = 1:num_segments,
    idx = (1:c.seg_size)+(i_p-1)*c.hop_size;       
    x = abs(fft(wav(idx).*c.w)/sum(c.w)*2).^2;
    
    P(:,i_p) = x(1:end/2+1);%工程实际中经常只用单边功率谱
   
end  

c.num_filt = 36; %% Mel频带数
f = linspace(0,c.fs/2,c.seg_size/2+1);%初始平均划分f
%%
mel = log(1+f/700)*1127.01048; %1127.01048=2595/log10 ,Matlab中log=ln

%%
mel_idx = linspace(0,mel(end),c.num_filt+2);%初始平均划分mel（38个点）
f_idx = zeros(c.num_filt+2,1);
for i=1:c.num_filt+2,
%% f_idx(i)存的是mel中与mel_idx(i)最近的元素的地址
   [tmp, f_idx(i)] = min(abs(mel - mel_idx(i)));%近似的平均划分       
end
%%
freqs = f(f_idx);
h = 2./(freqs(3:c.num_filt+2)-freqs(1:c.num_filt));%%三角的高度
c.mel_filter = zeros(c.num_filt,c.seg_size/2+1);

for i=1:c.num_filt,
  c.mel_filter(i,:) =(f > freqs(i) & f <= freqs(i+1)).* ...
                h(i).*(f-freqs(i))/(freqs(i+1)-freqs(i)) + ...
                (f > freqs(i+1) & f < freqs(i+2)).* ...
                h(i).*(freqs(i+2)-f)/(freqs(i+2)-freqs(i+1));
            M = zeros(c.num_filt,num_segments);  %初始化
end
for i_m = 1:num_segments,
    M(:,i_m) = c.mel_filter*P(:,i_m);% 通过三角滤波器
end

% 做对数变换
%M(M<1)=1;
M = 10*log(M);
% min(P(:))
% min(M(:))
% return
%DCT函数
c.DCT = 1/sqrt(c.num_filt/2) * ...
cos((0:num_ceps_coeffs-1)'*(0.5:c.num_filt)*pi/c.num_filt);
c.DCT(1,:) = c.DCT(1,:)*sqrt(2)/2;
c.DCT;
%%离散余弦变换
mfcc= c.DCT * M;
ncentres = 16;% 高斯分量个数
input_dim = 16; %特征维数
Z=mfcc';

% gmm(Z,16);

% % 设置混合模型
% mix = gmm(input_dim, ncentres, 'diag');
% % 特征数据输入
% siz=600;
% features = zeros(siz,input_dim);
%  for k=1:siz
%    for j=1:input_dim
%         features(k,j)=data.feat.mfcc(i_files,j,k);
%     end
% end
% % 初始化模型参数
% mix = gmminit(mix, features, options);
% options(14) = 20;% 迭代次数.
% [mix, options, errlog]=gmmem(mix, features, options);
% Gmmdata(i_files)=mix;

