%program to generate an exponentially damped sinusoidal wave


close all;
clear all;


A = input('Enter the amplitude of the sinusoidal wave A = '); %.15
f = input('Enter the frequency of the sinusoidal wave F = '); %100
phi = input('Enter the phase angle of sinusoid phi = '); %0
a = input('Enter the attenuation factor a = '); %90
i=1;
ecgptav = [];
for i = 1:30
%t = 0:.001:1;
t = 0:.001:.040;
ff = f+randn(size(t))*(f/10); 
w = ff*2*pi;
%y = abs(A*sin(w*t + phi).*exp(-a*t));
y = abs(A*sin(w + phi).*exp(-a*t));
y =[0 y 0];
% figure(1);
% plot(t,y);
%% Espectro de potencia
% Y = fft(y,size(y,2));
% Pyy = Y.*conj(Y)/size(y,2);
% F = 1000/size(y,2)*(0:42);
% figure(2)
% plot(F(1:23),Pyy(1:23))
% title('Power spectral density')
% xlabel('Frequency (Hz)')


%% Gernado o sinal ECG + PTAV
% ptav é adicionado nos pontos 211 até 251.
Y = y'; 
s = xlsread('sinalmedio');
z1 = zeros(209,1);
z2 = zeros(149,1);
ptav = [z1; Y; z2];
s1 = s + ptav;
% figure(3);
% plot(s1);
ecgptav(:,i) = s1;
i = i + 1;
end
%xlswrite('ECGPTAVG',ecgptav)

% YY = zeros(43,300);
% z = 0;
%  for i = 1:300;
%         ruido = 0.010;
%         yy = y + ruido*randn(size(y));
%         YY(:,i) = yy;
%         z = z + yy;
%  end
%  z = z/300;