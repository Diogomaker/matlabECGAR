%% Novo programa ECGAR
% Limpar workspace, fechar janelas e limpar janela de comandos
clear all; close all; clc;

% Importar o arquivo de ECG para an�lise de alta resolu��o
arquivo = uiimport;             %importa o arquivo para uma vari�vel de classe struct
vars = fieldnames(arquivo);     %mudan�a de classe para cell
sinal = arquivo.(vars{1});      %mudan�a de classe para matriz 

% Cada arquivo de ECG esta dividido da seguinte maneira: primeiro ponto dev X, segundo ponto dev Y e terceiro ponto dev Z e assim
%por diante em ciclos de tr�s.

%eliminar os tr�s primeiros pontos, pois o primeiro ponto de cada deriva��o � o valor
%de tens�o mais baixo da placa
sinal(1:3) = [];

% Alterando a escala de amplitude para milivolt:
% (conversor de 14 bits) dividir a tens�o por 8192
% (faixa de resolu��o 2^14 ou seja 16384 no total, metade na faixa positiva e metade na faixa negativa) 
% e multiplicar pelo ganho 5 mV 
sinal = (sinal/8192)*5;

x = []; % matriz da deriva��o X
y = []; % matriz da deriva��o Y
z = []; % matriz da deriva��o Z
n = length(sinal); %tamanho do vetor do sinal

% localiza��o dos pontos das deriva��es X,Y e Z:
a = 1; %contador
for i = 1:n/3
   x(i) = sinal(a);
   y(i) = sinal(a+1);
   z(i) = sinal(a+2);
   a = a + 3;
end
x = x';
y = y';
z = z';

% Filtragem passa-altas para remo��o de flutua��o de linha de base
fc = 0.5;       % frequ�ncia de corte
fn = 1000/2;    % frequ�ncia de nyquist = (frequ�ncia de amostragem)/2;
order = 2;      % segunda ordem
[b1, a1]=butter(order,(fc/fn),'high');
fx = filtfilt(b1,a1,x);  % utilizando a deriva��o X
%fx = filtfilt(b1,a1,y); % utilizando a deriva��o Y
%fx = filtfilt(b1,a1,z); % utilizando a deriva��o Z

%% Interpola��o para aumentar a frequ�ncia de amostragem
nf = 5;
fx = interp(fx,nf); % aumentando a frequencia de amostragem de 1kHz para 5kHz


%%  
% Pin�ando as amplitudes dos picos com uma dist�ncia de no m�nimo 2000 pontos entre eles
[maxwind,~] = findpeaks(fx, 'MINPEAKDISTANCE', 2000*nf );

% m�dia das amplitudes dos picos do sinal
windm = mean(maxwind);

%% Colocar na rotina op��es de entrada para definir o 'MINPEAKHEIGHT' e o 'MINPEAKDISTANCE' !!!!!!
%% Sele��o dos batimentos:
% Picos com valores m�nimo de 70% da m�dia dos picos, e dist�ncia entre os
% vizinhos de 700 pontos.
[pks,locs] = findpeaks(abs(fx), 'MINPEAKHEIGHT', windm*0.7, 'MINPEAKDISTANCE', 700*nf );

bat = fx(locs); % vetor de tens�es onde foram detectados batimentos
locs = locs';

% Gerar uma matriz com as informa��es :
% indice temporal para localiza��o do pico (primeira linha)
% valor RR entre o pico atual e o anterior (segunda linha)
MatInfo = [locs(2:end); diff(locs)]; % eliminado primeiro valor de batimento para igualar o tamanho dos vetores concatenados

% Matrix com os batimentos selecionados
batimentos = [];
for ii = 1 : length(MatInfo)-1
      batimentos(:,ii) = fx(locs(ii+1)-400*nf:locs(ii+1)+400*nf); % TAMANHO DOS BATIMENTOS SEPARADOS, a soma de 1 aos indices ii � para remo��o do primeiro ponto do vetor locs, pois o primeiro batimento possui tamanho comprometido
end
nbat1 = length(batimentos(1,:));

mediapontoref = mean((batimentos(400*nf+1,:))); 
% ponto m�ximo absoluto dos vetores batimentos � no indice 401

contador = [];
i = 0;
for ii = 1 : nbat1
    if (batimentos(400*nf+1,ii)) < mediapontoref*0.3 || (batimentos(400*nf+1,ii)) > mediapontoref*1.8
        i = i + 1;
        contador(i) = ii;     
     end
end
batimentos(:,contador) = []; % removendo os batimentos em que o indices de m�ximo absoluto fogem do padr�o da m�dia.
nbat2 = length(batimentos(1,:));

% Remo��o dos batimentos em que o picos fogem da m�dia 
mediapontomax = mean(max(abs(batimentos))); 
contador = [];
i = 0;
for ii = 1 : nbat2
    if max(abs(batimentos(:,ii))) < mediapontomax*0.3 || max(abs(batimentos(:,ii))) > mediapontomax*1.5
        i = i + 1;
        contador(i) = ii;     
     end
end
batimentos(:,contador) = []; % removendo os batimentos em que os picos de m�ximo absoluto fogem do padr�o da m�dia.
nbat3 = length(batimentos(1,:));
%% Calculo do ruido antes do alinhamento

% Voltar 220 pontos do pico do qrs e calcular o desvio padrao dos proximos
% 40 pontos, retirandoa a flutua��o da linha de base com a fun��o "detrend"
dp = [];
ruido_original = 0;
for i = 1: length(batimentos(1,:))
    [~,index] = max(abs(batimentos(:,i)));
    dp(i) = std(detrend(batimentos(index-220:index-180,i)));
    ruido_original = ruido_original + dp(i);
end
ruido_original = ruido_original/length(batimentos(1,:)); %Ruido do sinal original

%% Detector de n�vel

for i = 1:nbat3
    a1(i) = max(batimentos(:,i));
end
a2 = mean(a1)/2; % 50 por cento do valor da m�dia dos pontos m�ximos dos batimentos
ind = zeros(1,nbat3);
for i = 1:nbat3    
    [~, indmax] = max(abs(batimentos(:,i))) ; % indice do ponto m�ximo de cada batimento
    janela = abs(batimentos(indmax-50*nf:indmax+50*nf,i));
    t10up = find(janela >= a2, 1, 'first'); % primeiro indice com o valor mais pr�ximo ACIMA do nivel de 50% numa janela aberta em 100 pontos antes e depois do indmax
    t10down = t10up - 1; % primeiro indice com o valor mais pr�ximo ABAIXO do nivel de 50%
    if abs(janela(t10up) - a2) < abs(janela(t10down) - a2) % escolher qual indice � o mais pr�ximo de a2
        t10 = t10up;
    else
        t10 = t10down;
    end
    t11 = indmax - (50*nf+1-t10); % indice do nivel de 50% associada ao indmax que faz refer�ncia aos vetores batimentos
    ind(1,i) = t11;
    alinhado_dn(:,i) = batimentos((ind(1,i))-160:(ind(1,i))+180,i);     
end

%Sinal m�dio final (DN)
medio_dn = 0;
for i = 1:nbat3
   medio_dn = medio_dn + alinhado_dn(:,i);
end
medio_dn = medio_dn/nbat3;

%% Implementa��o do m�todo de Jan�, integral normalizada (NI - normalized integral)
% Um sinal integrado e noramlizado se assemelha morfologicamente � um
% degrau. A integral da diferen�a entre estes dois sinais integrados e
% normalizados � a diferen�a temporal entre eles.

% Alinhando os sinais de maneira que a integral da diferen�a entre o modelo
% e o sinal analisado seja nula.

%Integral Normalizada
DD = zeros(1,nbat3);
mod = batimentos(ind(1,1)-20*nf:ind(1,1)+60*nf,1); % modelo inicial
for i = 1:nbat3;
    modtri = mod.^2;
    modni = cumtrapz(modtri)/max(cumtrapz(modtri));
    tri = (batimentos(ind(1,i)-20*nf:ind(1,i)+60*nf,i)).^2;
    ni = cumtrapz(tri)/max(cumtrapz(tri));
    d = trapz(modni - ni);
    d = round(d);   % calcula d por arredondamento padr�o
    %d = ceil(d);   % calcula d por arredondamento para cima
    %d = floor(d);  % calcula d por arredondamento para baixo 
    DD(i) = d + ind(1,i);
    mod = mod.*(i-1);
    mod = mod + batimentos(ind(1,i)-20*nf+d:ind(1,i)+60*nf+d,i);
    mod = mod./i;
    alinhado_in(:,i) = batimentos(ind(1,i)-160+d:ind(1,i)+180+d,i);
end
    
%Sinal m�dio final (IN)
medio_in = 0;
for i = 1:nbat3
   medio_in = medio_in + alinhado_in(:,i);
end
medio_in = medio_in/nbat3;
    
%% Implemeta��o do m�todo de Integral Dupla

 % Integral Dupla
D = zeros(1,nbat3);
modelo = batimentos(ind(1,1)-20*nf:ind(1,1)+60*nf,1);
for i = 1:nbat3;
    lag = zeros(2,11);
    lag(1,:)=(-5:5);
    modeloid = cumtrapz((cumtrapz((modelo).^2)).^2);
    for l = -5:+5
        id = cumtrapz((cumtrapz((batimentos(ind(1,i)-20*nf+l:ind(1,i)+60*nf+l,i).^2))).^2);
        %id = batimentos(T1-20+l:T1+60+l,i);
        cc = corrcoef(id,modeloid);
        lag(2,l+6) = cc(1,2);
    end
    [c,k] = max(lag(2,:));
    d = lag(1,k);
    D(i) = d + ind(1,i);
    modelo = modelo.*(i-1);
    modelo = modelo + batimentos(ind(1,i)-20*nf+d:ind(1,i)+60*nf+d,i);
    modelo = modelo./i;
    alinhado_id(:,i) = batimentos(ind(1,i)-160+d:ind(1,i)+180+d,i);
end
%Sinal m�dio final (ID)
medio_id = 0;
for i = 1:nbat3
    medio_id = medio_id + alinhado_id(:,i);
end
medio_id = medio_id/nbat3;

figure
hold on
plot(medio_id,'b')
plot(medio_in,'r')
plot(medio_dn,'m')
legend('id','in','dn')
title('sinais medios')

figure
plot(batimentos(400*nf+1-170*nf:400*nf+1+170*nf+1,:))
title('batimentos')

%% pensar em um jeito para escolher qual o m�todo ser� usado, ou se os tres ser�o usados e depois escolher o resultado....
% arranjar um meio de salvar as informa��es finais em um arquivo; 
% quais informa��es salvar:
% ruidos, erro de alinhamento, potencia espectral
%
% utilizar sinais de individuos com as mesmas caracteristicas 
% no grupo dente e controle

%% �ndices de desempenho
% ruido final
[~,index] = max(medio_dn);
ruido_final_dn = std(detrend(medio_dn(index-71:index-41,1)));

[~,index] = max(medio_in);
ruido_final_in = std(detrend(medio_in(index-71:index-41,1)));

[~,index] = max(medio_id);
ruido_final_id = std(detrend(medio_id(index-71:index-41,1)));

% erro de alinhamento
errodn = std(ind(1,:));
erroin = std(DD);
erroid = std(D);


%% Tabela resposta

metodos = {'Sinal original';'Detector de n�vel'; 'Integral Normalizada'; 'Integral Dupla'};
ruidos = [ruido_original; ruido_final_dn; ruido_final_in; ruido_final_id];
erro_alinhamento = [0; errodn; erroin; erroid];
T = table(ruidos,erro_alinhamento,...
    'RowNames',metodos)