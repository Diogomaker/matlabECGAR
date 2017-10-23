%% Espectro de potência
%Calcular fft do sinal modelo e do sinal médio para obtenção da frequência de corte do filtro passa baixa equivalente e 
% comparar com a frequência de corte obtida pela formula “w0 = 125/desvio padrão (EA)”.
%dados utilizados:
%s(ids+18:ids+58,1)
%medio_dn(iddn+difi:iddn+diff,1))
%medio_in(idin+difi:idin+diff,1))
%medio_id(idid+difi:idid+diff,1))
y1 = input('Escolha o vetor a ser calculado o espectro = ');
%y1 = medio_in(idin+difi:idin+diff,1)'; %sinal a ser analisado 1, obs.: precisa ser vetor linha.

% fa = 1000; %frequencia de amostragem
% Y1 = fft(y1,size(y1,2));
% Pyy1 = Y1.*conj(Y1)/size(y1,2);
% F = fa/size(y1,2)*(0:42); 
% 
% figure
% plot(F(1:23),Pyy1(1:23))
% title('Densidade espectral de potência')
% xlabel('Frenquência (Hz)')
% 
% figure
% semilogy(F(1:23),Pyy1(1:23))
% xlabel('Frenquência (Hz)','fontsize',18)
% ylabel('Potência', 'fontsize',18)
% grid on

%% USANDO JANELA DE HANNING
fa = 1000;
L = length(y1);
La = round(L/2);
winvec = hann(L); 
Y1 = fft(y1'.*winvec);
Pyy1 = Y1.*conj(Y1)/L;
F = fa/L*(0:La); 
plot(F(1:La/8),Pyy1(1:La/8))
grid on
%figure
%semilogy(F(1:La),Pyy1(1:La))
xlabel('Frenquência (Hz)','fontsize',18)
ylabel('Potência', 'fontsize',18)
