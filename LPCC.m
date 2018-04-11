% close all
% clear all
% clc
% MM=307244;
% I = audioread('work50-1-all 00_00_00-00_00_30.wav',[1 MM]);%读入原始语音
function Z=LPCC(wav)
I=wav(:,1);
% subplot(3,1,1),plot(I);
% title('原始语音波形')
%对指定帧位置进行加窗处理
Q = I';
N = 256; % 窗长
Hamm = hamming(N); % 加窗
frame = 30;%需要处理的帧位置
M = Q(((frame - 1) * (N / 2) + 1):((frame - 1) * (N / 2) + N));
Frame = M .* Hamm';%加窗后的语音帧
[B,F,T] = specgram(I,N,N/2,N);
[m,n] = size(B);
for i = 1:m
FTframe1(i) = B(i,frame);
end
%%%%%
P =12; % 预测器阶数

ai = lpc(Frame,P); % 计算lpc系数
LP = filter([0 -ai(2:end)],1,Frame); % 建立语音帧的正则方程
FFTlp = fft(LP);
E = Frame - LP; % 预测误差
% subplot(3,1,2),plot(1:N,Frame,1:N,LP,'-r');grid;
% title('原始语音和预测语音波形')
% subplot(3,1,3),plot(E);grid;
% title('预测误差');

fLength(1 : 2 * N) = [M,zeros(1,N)];
Xm = fft(fLength,2 * N);
X = Xm .* conj(Xm);
Y = fft(X , 2 * N);
Rk = Y(1 : N);
PART = sum(ai(2 : P + 1) .* Rk(1 : P));
G = sqrt(sum(Frame.^2) - PART);
A = (FTframe1 - FFTlp(1 : length(F'))) ./ FTframe1 ;
% subplot(2,1,1),plot(F',20*log(abs(FTframe1)),F',(20*log(abs(1 ./ A))),'-r');grid;
% xlabel('频率/dB');ylabel('幅度');
% title('短时谱');
% subplot(2,1,2),plot(F',(20*log(abs(G ./ A))));grid;
% xlabel('频率/dB');ylabel('幅度');
% 
% title('LPC谱');



%求出预测误差的倒谱

pitch = fftshift(rceps(E));

M_pitch = fftshift(rceps(Frame));
Z=M_pitch;

% subplot(2,1,1),plot(M_pitch);grid;
% 
% xlabel('语音帧');ylabel('/dB');
% 
% title('原始语音帧倒谱');
% 
% subplot(2,1,2),plot(pitch);grid;
% 
% xlabel('语音帧');ylabel('/dB');
% 
% title('预测误差倒谱');
% 
% 
% 
% %画出语谱图
% 
% ai1 = lpc(I,P); % 计算原始语音lpc系数
% 
% LP1 = filter([0 -ai(2:end)],1,I); % 建立原始语音的正则方程
% 
% subplot(2,1,1);
% 
% specgram(I,N,N/2,N);
% 
% title('原始语音语谱图');
% 
% subplot(2,1,2);
% 
% specgram(LP1,N,N/2,N);
% 
% title('预测语音语谱图');
