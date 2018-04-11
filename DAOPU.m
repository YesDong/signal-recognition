clc
clear all
close all
ts = 0.001;
Fs = 1/ts;
t = 0:ts:1;
L = length(t);
f0 = 100;
fm = 5;
A = 1;
% Ô­Ê¼ÐÅºÅ
y = A*(1+cos(2*pi*fm*t)+cos(2*pi*2*fm*t)+cos(2*pi*3*fm*t)+cos(2*pi*4*fm*t)+cos(2*pi*5*fm*t)).*sin(2*pi*f0*t);
figure
plot(t, y)
% ÆµÆ×
NFFT = 2^nextpow2(L);
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure
plot(f,2*abs(Y(1:NFFT/2+1)))
% µ¹ÆµÆ×
c = rceps(y);
figure
plot(t,c)