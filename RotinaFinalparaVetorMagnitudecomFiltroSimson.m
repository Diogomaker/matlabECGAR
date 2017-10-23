%% ECGAR para obter vetor magnitude com filtro bidirecional de SIMSON, 1981
% Limpar workspace, fechar janelas e limpar janela de comandos
clear all; close all; clc;
%Frequencia de amostragem = 1000Hz

% Importar o arquivo de ECG para análise de alta resolução
arquivo = uiimport;             %importa o arquivo para uma variável de classe struct
vars = fieldnames(arquivo);     %mudança de classe para cell
sinal = arquivo.(vars{1});      %mudança de classe para matriz 

% Cada arquivo de ECG esta dividido da seguinte maneira: primeiro ponto dev X, segundo ponto dev Y e terceiro ponto dev Z e assim
%por diante em ciclos de três.

%eliminar os três primeiros pontos, pois o primeiro ponto de cada derivação é o valor
%de tensão mais baixo da placa
sinal(1:3) = [];

% Alterando a escala de amplitude para milivolt:
% (conversor de 14 bits) dividir a tensâo por 8192
% (faixa de resolução 2^14 ou seja 16384 no total, metade na faixa positiva e metade na faixa negativa) 
% e multiplicar pelo ganho 5 mV 
sinal = (sinal/8192)*5;

x = []; % matriz da derivação X
y = []; % matriz da derivação Y
z = []; % matriz da derivação Z
n = length(sinal); %tamanho do vetor do sinal

% localização dos pontos das derivações X,Y e Z:
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


% Filtragem passa-altas para remoção de flutuação de linha de base
fc = 0.5;       % frequência de corte
fn = 1000/2;    % frequência de nyquist = (frequência de amostragem)/2;
order = 2;      % segunda ordem
[b1, a1]=butter(order,(fc/fn),'high');
fxx = filtfilt(b1,a1,x);  % utilizando a derivação X
fy = filtfilt(b1,a1,y); % utilizando a derivação Y
fz = filtfilt(b1,a1,z); % utilizando a derivação Z

%% !!!CALCULAR O VALOR MEDIO ALINHADOS DAS DERIVAÇÕES X Y E Z E DEPOIS FAZER A FILTRAGEM PASSABANDA E PEGAR O VALOR ABSOLUTO!!!

f = [fxx fy fz];
%%
Medio_dn = [];
Medio_in = [];
Medio_id = [];

for iii = 1:3
    fx = f(:,iii);
    % Pinçando as amplitudes dos picos com uma distância de no mínimo 2000 pontos entre eles
    [maxwind,~] = findpeaks(fx, 'MINPEAKDISTANCE', 2000 );

    % média das amplitudes dos picos do sinal
    windm = mean(maxwind);

    %% Colocar na rotina opções de entrada para definir o 'MINPEAKHEIGHT' e o 'MINPEAKDISTANCE' !!!!!!
    %% Seleção dos batimentos:
    % Picos com valores mínimo de 70% da média dos picos, e distância entre os
    % vizinhos de 700 pontos.
    [pks,locs] = findpeaks(abs(fx), 'MINPEAKHEIGHT', windm*0.7, 'MINPEAKDISTANCE', 700 );

    bat = fx(locs); % vetor de tensões onde foram detectados batimentos
    locs = locs';

    % Gerar uma matriz com as informações :
    % indice temporal para localização do pico (primeira linha)
    % valor RR entre o pico atual e o anterior (segunda linha)
    MatInfo = [locs(2:end); diff(locs)]; % eliminado primeiro valor de batimento para igualar o tamanho dos vetores concatenados

    % Matrix com os batimentos selecionados
    batimentos = [];
    for ii = 1 : length(MatInfo)-1
          batimentos(:,ii) = fx(locs(ii+1)-400:locs(ii+1)+400); % TAMANHO DOS BATIMENTOS SEPARADOS, a soma de 1 aos indices ii é para remoção do primeiro ponto do vetor locs, pois o primeiro batimento possui tamanho comprometido
    end
    nbat1 = length(batimentos(1,:));

    mediapontoref = mean(abs(batimentos(401,:))); 
    % ponto máximo absoluto dos vetores batimentos é no indice 401

    contador = [];
    i = 0;
    for ii = 1 : nbat1
        if abs(batimentos(401,ii)) < abs(mediapontoref)*0.5 || abs(batimentos(401,ii)) > abs(mediapontoref)*1.8
            i = i + 1;
            contador(i) = ii;     
         end
    end
    batimentos(:,contador) = []; % removendo os batimentos em que o indices de máximo absoluto fogem do padrão da média.
    nbat2 = length(batimentos(1,:));

    % Remoção dos batimentos em que o picos fogem da média 
    mediapontomax = mean(max(abs(batimentos))); 
    contador = [];
    i = 0;
    for ii = 1 : nbat2
        if max(abs(batimentos(:,ii))) < mediapontomax*0.3 || max(abs(batimentos(:,ii))) > mediapontomax*1.5
            i = i + 1;
            contador(i) = ii;     
         end
    end
    % removendo os batimentos em que os picos de máximo absoluto fogem do padrão da média.
    batimentos(:,contador) = []; 
    nbat3 = length(batimentos(1,:));
    
    %% Detector de nível parte 1 (sem aplicar corralação)

    for i = 1:nbat3
        a1(i) = max(abs(batimentos(:,i)));
    end
    a2 = mean(a1)/2; % 50 por cento do valor da média dos pontos máximos dos batimentos
    ind = zeros(1,nbat3);
    for i = 1:nbat3    
%         [~, indmax] = max(abs(batimentos(:,i))) ; % indice do ponto máximo de cada batimento
        indmax = 401;
        janela = abs(batimentos(indmax-70:indmax+70,i));
        t10up = find(janela >= a2, 1, 'first'); % primeiro indice com o valor mais próximo ACIMA do nivel de 50% numa janela aberta em 100 pontos antes e depois do indmax
        if t10up == 1; % passar para a proxima iteração no caso do t10up ser 1 
            t10 = t10up;
            t11 = indmax - (71-t10); % indice do nivel de 50% associada ao indmax que faz referência aos vetores batimentos
            ind(1,i) = t11;
            alinhado_dn(:,i) = batimentos((ind(1,i))-160:(ind(1,i))+180,i);
            continue
        end
        t10down = t10up - 1; % primeiro indice com o valor mais próximo ABAIXO do nivel de 50%
        if abs(janela(t10up) - a2) < abs(janela(t10down) - a2) % escolher qual indice é o mais próximo de a2
            t10 = t10up;
        else
            t10 = t10down;
        end
        t11 = indmax - (71-t10); % indice do nivel de 50% associada ao indmax que faz referência aos vetores batimentos
        ind(1,i) = t11;
        alinhado_dn(:,i) = batimentos((ind(1,i))-160:(ind(1,i))+180,i);     
    end

    %Sinal médio final (DN)
    medio_dn = 0;
    for i = 1:nbat3
       medio_dn = medio_dn + alinhado_dn(:,i);
    end
    medio_dn = medio_dn/nbat3;
    %% Removendo os batimentos com baixa correlação
    for i = 1:nbat3
        co(i) = corr(batimentos((ind(1,i))-160:(ind(1,i))+180,i),medio_dn);
    end
    refco = mean(co);
    contador = [];
    i = 0;
    for ii = 1 : nbat3
        if co(ii) <  refco
            i = i + 1;
            contador(i) = ii;     
         end
    end
    batimentos(:,contador) = []; % removendo os batimentos em que os picos de máximo absoluto fogem do padrão da média.
    nbat4 = length(batimentos(1,:));
    
    %% Filtragem dos batimentos para indexação dos pontos de referência
    % os sinais filtrados serão utilizados para deteção dos pontos de
    % alinhamento, porém os sinais alinhados serão os sinais sem filtragem
    % passa-baixas
    fc = 40;       % frequência de corte
    fn = 1000/2;    % frequência de nyquist = (frequência de amostragem)/2;
    order = 2;      % segunda ordem
    [b1, a1]=butter(order,(fc/fn),'low');
    batimentosf = filtfilt(b1,a1,batimentos);
    %% Detector de nível parte 2
    for i = 1:nbat4
        a1(i) = max(abs(batimentosf(:,i)));
    end
    a2 = mean(a1)/2; % 50 por cento do valor da média dos pontos máximos dos batimentos
    ind = zeros(1,nbat4);
    for i = 1:nbat4   
        [~, indmax] = max(abs(batimentosf(:,i))) ; % indice do ponto máximo de cada batimento
        janela = abs(batimentosf(indmax-50:indmax+50,i));
        t10up = find(janela >= a2, 1, 'first'); % primeiro indice com o valor mais próximo ACIMA do nivel de 50% numa janela aberta em 100 pontos antes e depois do indmax
        if t10up == 1; % passar para a proxima iteração no caso do t10up ser 1 
            t10 = t10up;
            t11 = indmax - (71-t10); % indice do nivel de 50% associada ao indmax que faz referência aos vetores batimentos
            ind(1,i) = t11;
            alinhado_dn(:,i) = batimentos((ind(1,i))-160:(ind(1,i))+180,i);
            continue
        end
        t10down = t10up - 1; % primeiro indice com o valor mais próximo ABAIXO do nivel de 50%
        if abs(janela(t10up) - a2) < abs(janela(t10down) - a2) % escolher qual indice é o mais próximo de a2
            t10 = t10up;
        else
            t10 = t10down;
        end
        t11 = indmax - (51-t10); % indice do nivel de 50% associada ao indmax que faz referência aos vetores batimentos
        ind(1,i) = t11;
        alinhado_dn(:,i) = batimentos((ind(1,i))-160:(ind(1,i))+180,i);     
    end

    %Sinal médio final (DN)
    medio_dn = 0;
    for i = 1:nbat4
       medio_dn = medio_dn + alinhado_dn(:,i);
    end
    medio_dn = medio_dn/nbat4;

    %% Implementação do método de Jané, integral normalizada (NI - normalized integral)
    % Um sinal integrado e noramlizado se assemelha morfologicamente à um
    % degrau. A integral da diferença entre estes dois sinais integrados e
    % normalizados é a diferença temporal entre eles.

    % Alinhando os sinais de maneira que a integral da diferença entre o modelo
    % e o sinal analisado seja nula.

    %Integral Normalizada
    DD = zeros(1,nbat4);
    mod = batimentosf(ind(1,1)-20:ind(1,1)+60,1); % modelo inicial
    for i = 1:nbat4;
        modtri = mod.^2;
        modni = cumtrapz(modtri)/max(cumtrapz(modtri));
        tri = (batimentosf(ind(1,i)-20:ind(1,i)+60,i)).^2;
        ni = cumtrapz(tri)/max(cumtrapz(tri));
        d = trapz(modni - ni);
        d = round(d);   % calcula d por arredondamento padrão
        %d = ceil(d);   % calcula d por arredondamento para cima
        %d = floor(d);  % calcula d por arredondamento para baixo 
        DD(i) = d + ind(1,i);
        mod = mod.*(i-1);
        mod = mod + batimentosf(ind(1,i)-20+d:ind(1,i)+60+d,i);
        mod = mod./i;
        alinhado_in(:,i) = batimentos(ind(1,i)-160+d:ind(1,i)+180+d,i);
    end

    %Sinal médio final (IN)
    medio_in = 0;
    for i = 1:nbat4
       medio_in = medio_in + alinhado_in(:,i);
    end
    medio_in = medio_in/nbat4;

    %% Implemetação do método de Integral Dupla

     % Integral Dupla
    D = zeros(1,nbat4);
    modelo = batimentosf(ind(1,1)-20:ind(1,1)+60,1);
    for i = 1:nbat4;
        lag = zeros(2,11);
        lag(1,:)=(-5:5);
        modeloid = cumtrapz((cumtrapz((modelo).^2)).^2);
        for l = -5:+5
            id = cumtrapz((cumtrapz((batimentosf(ind(1,i)-20+l:ind(1,i)+60+l,i).^2))).^2);
            %id = batimentosf(T1-20+l:T1+60+l,i);
            cc = corrcoef(id,modeloid);
            lag(2,l+6) = cc(1,2);
        end
        [c,k] = max(lag(2,:));
        d = lag(1,k);
        D(i) = d + ind(1,i);
        modelo = modelo.*(i-1);
        modelo = modelo + batimentosf(ind(1,i)-20+d:ind(1,i)+60+d,i);
        modelo = modelo./i;
        alinhado_id(:,i) = batimentos(ind(1,i)-160+d:ind(1,i)+180+d,i);
    end
    %Sinal médio final (ID)
    medio_id = 0;
    for i = 1:nbat4
        medio_id = medio_id + alinhado_id(:,i);
    end
    medio_id = medio_id/nbat4;

%     figure
%     hold on
%     plot(medio_id,'b')
%     plot(medio_in,'r')
%     plot(medio_dn,'m')
%     legend('id','in','dn')
%     title('sinais medios')
% 
%     figure
%     plot(batimentos(401-170:401+171,:))
%     title('batimentos')

    %% pensar em um jeito para escolher qual o método será usado, ou se os tres serão usados e depois escolher o resultado....
    % arranjar um meio de salvar as informações finais em um arquivo; 
    % quais informações salvar:
    % ruidos, erro de alinhamento, potencia espectral
    %
    % utilizar sinais de individuos com as mesmas caracteristicas 
    % no grupo dente e controle

    %% Ìndices de desempenho
    % ruido final
%     [~,index] = max(medio_dn);
%     ruido_final_dn = std(detrend(medio_dn(index-71:index-41,1)));
% 
%     [~,index] = max(medio_in);
%     ruido_final_in = std(detrend(medio_in(index-71:index-41,1)));
% 
%     [~,index] = max(medio_id);
%     ruido_final_id = std(detrend(medio_id(index-71:index-41,1)));

    % erro de alinhamento
%     errodn = std(ind(1,:));
%     erroin = std(DD);
%     erroid = std(D);

    if iii == 1
        X = [medio_dn medio_in medio_id];
    elseif iii == 2
        Y = [medio_dn medio_in medio_id];
    elseif iii == 3
        Z = [medio_dn medio_in medio_id];
    end
end
%% Grampeamento para evitar ringing nos sinais filtrados 
% X
aX = X(1,:);
bX = X(end,:);
L = length(X);

xx = (1:1:L);
Xgra = [];
xbase = [1 L];
ybase = [aX(1) bX(1)];

a = polyfit(xbase,ybase,1);
reta = a(1)*xx+a(2);
for i =1:3
    Xgra(:,i) = X(:,i) - reta';
end

% Y
aY = Y(1,:);
bY = Y(end,:);
L = length(Y);

yy = (1:1:L);
Ygra = [];
xbase = [1 L];
ybase = [aY(1) bY(1)];

a = polyfit(xbase,ybase,1);
reta = a(1)*yy+a(2);
for i =1:3
    Ygra(:,i) = Y(:,i) - reta';
end

%Z
aZ = Z(1,:);
bZ = Z(end,:);
L = length(Z);

zz = (1:1:L);
Zgra = [];
xbase = [1 L];
ybase = [aZ(1) bZ(1)];

a = polyfit(xbase,ybase,1);
reta = a(1)*zz+a(2);
for i =1:3
    Zgra(:,i) = Z(:,i) - reta';
end


%% Filtragem para geração do VM
% os sinais filtrados serão utilizados para deteção dos pontos de
% alinhamento, porém os sinais alinhados serão os sinais sem filtragem
% passa-baixas.
%     De uma maneira geral, na análise do ECGAR no Domínio do Ttempo,
% utiliza-se o filtro digital Butterworth de 4 pólos passa-faixas, com frequência de corte
% passa-altas de 25, 40 ou 80 Hz e passa-baixas em 250 Hz (SIMSON, 1981). Em 1981,
% Simson desenvolveu o filtro que se tornou padrão para análise dos PTAV, o qual é
% aplicado de maneira bidirecional (dos extremos do sinal medio para um ponto médio, no
% interior do complexo QRS) para evitar a distorção de fase nos componentes de
% frequência do sinal próximo aos valores de corte.
 
fcutlow = 40;       
fcuthigh = 250;
fn = 1000/2;    % frequência de nyquist = (frequência de amostragem)/2;
order = 4;      % segunda ordem
[b1, a1] = butter(order,[fcutlow,fcuthigh]/(fn), 'bandpass');

FX1 = filter(b1,a1,Xgra(1:170,:)); % Filtragem até o meio do QRS
Xinv = flipud(Xgra); % batimentos invertidos
FX2inv = filter(b1,a1,Xinv(1:171,:));
FX2 = flipud(FX2inv); 
FX = [FX1; FX2]; 

FY1 = filter(b1,a1,Ygra(1:170,:)); % Filtragem até o meio do QRS
Yinv = flipud(Ygra); % batimentos invertidos
FY2inv = filter(b1,a1,Yinv(1:171,:));
FY2 = flipud(FY2inv); 
FY = [FY1; FY2];

FZ1 = filter(b1,a1,Zgra(1:170,:)); % Filtragem até o meio do QRS
Zinv = flipud(Zgra); % batimentos invertidos
FZ2inv = filter(b1,a1,Zinv(1:171,:));
FZ2 = flipud(FZ2inv); 
FZ = [FZ1; FZ2];

AFX = abs(FX);
AFY = abs(FY);
AFZ = abs(FZ);
vm_dn = sqrt(AFX(:,1).^2+AFY(:,1).^2+AFZ(:,1).^2);
vm_in = sqrt(AFX(:,2).^2+AFY(:,2).^2+AFZ(:,2).^2);
vm_id = sqrt(AFX(:,3).^2+AFY(:,3).^2+AFZ(:,3).^2);

%% Rotina para calcular dQRS, LAS40 e RMS40; 
% Usar depois de "RotinaFinalparaVetorMagnitude.m"
% dQRS: duração da ativação ventricular
% LAS40: duração do segmento terminal do complexo QRS com amplitude abaixo de 40 uV
% RMS40: média quadrática da amplitude nos 40 ms terminais do complexo QRS
%% DN
Vm = vm_dn;
N = length (AFX);

[magnitude, indmag] = max(Vm);        % Valor absoluto do ponto maximo
% calcular o rms de uma janela movel de 40 pontos até o ponto de max
% amplitude

% Medir a variancia em janelas de 40 pontos do começo até o final do
% batimento. Com isso a partir do ponto de menor variancia ir calculando a
% variancia em janelas de 5 pontos até achar tres janelas consecutivas com
% a variancia 5x maior do que a variancia minima


% Final QRS
var40 = [];
for i = indmag:length(Vm)-40
    jan40 = Vm(i:i+40);
    var40(i) = var(jan40);
end
% calcular a variancia minima das janelas de 40 pontos e seu indice para servir
% de parametro para a busca o inicio do qrs
[varmin,indmin] = min(var40(indmag:length(Vm)-40));
indmin = indmag - 1 + indmin;
var5 = [];
for i = indmin: -1:indmag
    jan5 = Vm(i-5:i);
    var5(i) = var(jan5);
end
achados = find(var5>5*varmin);
tam = length(achados);
i = tam; j = 0;
while j == 0
    if achados(i) == achados(i-1)+1 && achados(i-1)+1 == achados(i-2)+2
        qrs_fim = achados(i);
        j = 1;
    else
        i = i-1;
    end
end

% Inicio QRS
var5 = [];
for i = 75:indmag
    jan5 = Vm(i:i+5);
    var5(i) = var(jan5);
end
[varmin,indmin] = min(var5(75:indmag));
indmin = 75 - 1 + indmin;
[~,imini] = min(Vm(75:indmag));
imini = 75 - 1 + imini;

imin = indmin;
% if imini >= indmin
%     imin = imini;
% else
%     imin = indmin;
% end

achados = find(var5(imin:end)>5*varmin);
achados = achados + imin-1;
tam = length(achados);

i = 1; j = 0;
while j == 0
    if achados(i) == achados(i+1)-1 && achados(i+1)-1 == achados(i+2)-2
        qrs_inicio = achados(i);
        j = 1;
    else
        i = i+1;
    end
end

qrs_fim_dn = qrs_fim;
qrs_inicio_dn = qrs_inicio;
dQRS_dn = qrs_fim_dn - qrs_inicio_dn;

% RMS40 (valor RMS dos 40 ms terminais da ativaçao ventricular)
u = 1;
for i = qrs_fim_dn - 40 : qrs_fim_dn;
    QRS_40_dn(u) = Vm(i);
    u = u + 1;
end;

RMS40_dn = (sqrt(sum(QRS_40_dn.^2)/length(QRS_40_dn)))*1000;



% LAS40 (duraçao dos potenciais com amplitude inferior a 40 microV na regiao terminal do QRS)
 LAS40_dn = 0;
i = qrs_fim_dn;
while (Vm(i) <= 0.04);
    LAS40_dn = LAS40_dn + 1;
    i = i - 1;
end

LAS40_dn = LAS40_dn - 1 ;


%% IN
Vm = vm_in;
N = length (AFX);

[magnitude,indmag] = max(Vm);        % Valor absoluto do ponto maximo
% Final QRS
var40 = [];
for i = indmag:length(Vm)-40
    jan40 = Vm(i:i+40);
    var40(i) = var(jan40);
end
% calcular a variancia minima das janelas de 40 pontos e seu indice para servir
% de parametro para a busca o inicio do qrs
[varmin,indmin] = min(var40(indmag:length(Vm)-40));
indmin = indmag - 1 + indmin;
var5 = [];
for i = indmin: -1:indmag
    jan5 = Vm(i-5:i);
    var5(i) = var(jan5);
end
achados = find(var5>5*varmin);
tam = length(achados);
i = tam; j = 0;
while j == 0
    if achados(i) == achados(i-1)+1 && achados(i-1)+1 == achados(i-2)+2
        qrs_fim = achados(i);
        j = 1;
    else
        i = i-1;
    end
end

% Inicio QRS
var5 = [];
for i = 75:indmag
    jan5 = Vm(i:i+5);
    var5(i) = var(jan5);
end
[varmin,indmin] = min(var5(75:indmag));
indmin = 75 - 1 + indmin;
[~,imini] = min(Vm(75:indmag));
imini = 75 - 1 + imini;

imin = indmin;
% if imini >= indmin
%     imin = imini;
% else
%     imin = indmin;
% end

achados = find(var5(imin:end)>5*varmin);
achados = achados + imin-1;
tam = length(achados);

i = 1; j = 0;
while j == 0
    if achados(i) == achados(i+1)-1 && achados(i+1)-1 == achados(i+2)-2
        qrs_inicio = achados(i);
        j = 1;
    else
        i = i+1;
    end
end


qrs_fim_in = qrs_fim;
qrs_inicio_in = qrs_inicio;
dQRS_in = qrs_fim_in - qrs_inicio_in;

% RMS40 (valor RMS dos 40 ms terminais da ativaçao ventricular)
u = 1;
for i = qrs_fim_in - 40 : qrs_fim_in;
    QRS_40_in(u) = Vm(i);
    u = u + 1;
end;

RMS40_in = (sqrt(sum(QRS_40_in.^2)/length(QRS_40_in)))*1000;



% LAS40 (duraçao dos potenciais com amplitude inferior a 40 microV na regiao terminal do QRS)
 LAS40_in = 0;
i = qrs_fim_in;
while (Vm(i) <= 0.04);
    LAS40_in = LAS40_in + 1;
    i = i - 1;
end

LAS40_in = LAS40_in - 1 ;


%% ID
Vm = vm_id;
N = length (AFX);

[magnitude,indmag] = max(Vm);        % Valor absoluto do ponto maximo
% Final QRS
var40 = [];
for i = indmag:length(Vm)-40
    jan40 = Vm(i:i+40);
    var40(i) = var(jan40);
end
% calcular a variancia minima das janelas de 40 pontos e seu indice para servir
% de parametro para a busca o inicio do qrs
[varmin,indmin] = min(var40(indmag:length(Vm)-40));
indmin = indmag - 1 + indmin;
var5 = [];
for i = indmin: -1:indmag
    jan5 = Vm(i-5:i);
    var5(i) = var(jan5);
end
achados = find(var5>5*varmin);
tam = length(achados);
i = tam; j = 0;
while j == 0
    if achados(i) == achados(i-1)+1 && achados(i-1)+1 == achados(i-2)+2 && achados(i-2)+2 == achados(i-3)+3
        qrs_fim = achados(i);
        j = 1;
    else
        i = i-1;
    end
end

% Inicio QRS
var5 = [];
for i = 75:indmag
    jan5 = Vm(i:i+5);
    var5(i) = var(jan5);
end
[varmin,indmin] = min(var5(75:indmag));
indmin = 75 - 1 + indmin;
[~,imini] = min(Vm(75:indmag));
imini = 75 - 1 + imini;

imin = indmin;
% if imini >= indmin
%     imin = imini;
% else
%     imin = indmin;
% end

achados = find(var5(imin:end)>5*varmin);
achados = achados + imin-1;
tam = length(achados);

i = 1; j = 0;
while j == 0
    if achados(i) == achados(i+1)-1 && achados(i+1)-1 == achados(i+2)-2
        qrs_inicio = achados(i);
        j = 1;
    else
        i = i+1;
    end
end


qrs_fim_id = qrs_fim;
qrs_inicio_id = qrs_inicio;
dQRS_id = qrs_fim_id - qrs_inicio_id;

% RMS40 (valor RMS dos 40 ms terminais da ativaçao ventricular)
u = 1;
for i = qrs_fim_id - 40 : qrs_fim_id;
    QRS_40_id(u) = Vm(i);
    u = u + 1;
end;

RMS40_id = (sqrt(sum(QRS_40_id.^2)/length(QRS_40_id)))*1000;



% LAS40 (duraçao dos potenciais com amplitude inferior a 40 microV na regiao terminal do QRS)
 LAS40_id = 0;
i = qrs_fim_id;
while (Vm(i) <= 0.04);
    LAS40_id = LAS40_id + 1;
    i = i - 1;
end

LAS40_id = LAS40_id - 1 ;


% Tabela de resultados
cab = {'dQRS (ms)';'RMS40 (uV)';'LAS40 (ms)'};
dn = [dQRS_dn; RMS40_dn; LAS40_dn];
in = [dQRS_in; RMS40_in; LAS40_in];
id = [dQRS_id; RMS40_id; LAS40_id];
Tabela = table(dn,in,id,...
    'RowNames',cab)