% Rotina para calcular dQRS, LAS40 e RMS40; 
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