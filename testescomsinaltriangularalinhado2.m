%% Alinhando os sinais com DN, ID , NI , DN+ID e DN+NI 

clear all; clc; close all;

%filtragem passa-baixas
FC = 60/500;
[B,A]= butter(2,FC);
for ii = 1:100
    batimentos = zeros(300,400);
    batimentosf = zeros(300,400);
    for i = 1:400;
        T = 1*(1/12.5);
        Fs = 1000;
        dt = 1/Fs;
        t = 0:dt:T-dt;
        x = sawtooth(2*pi*12.5*t,0.5);
        x = (x+1)/2;
        z = zeros(1,109);
        x = [z x z 0 0];
        s = x';
        ruido = 0.090;
        r = ruido*randn(size(x));
        x = x + r;
        x = x';
        batimentosf(:,i) = x;
        y = filtfilt(B,A,x);
        batimentos(:,i) = y;
    end
    [~,pico] = max(s);
    A2 = max(s)/2; 
    T1 = find(s >= A2 ,1,'first');
    [~,T2] = max(s); 
    
    for i = 1:400
       [~,pc] = max(batimentos(:,i));
       picos(i) = pc;
    end
    
    % DETECTOR DE NÍVEL

    %     tempos = zeros(3,400);
    %     for i = 1:400
    %         a2 = max(batimentos(:,i))/2; 
    %         t1 = find(batimentos(:,i)>=a2,1,'first');
    %         t2 = find(batimentos(:,i)>=a2,1,'last');
    %         tempos(1,i) = t1;    
    %         tempos(2,i) = t2;   
    %         tempos(3,i) = (t1+t2)/2;        
    %     end
    
    tempos = zeros(3,400);
    for i = 1:400
        a1(i) = max(batimentos(:,i));
    end
    
    for i = 1:400
        [a indmax] = max(batimentos(:,i)); 
        %a2 = a/2; %valor procurado 
        a2 = mean(a1)/2;
        tmp1 = abs(batimentos(1:indmax,i)-a2); %valores da diferença até ponto maximo
        tmp = [zeros(indmax-1,1); batimentos(indmax:end,i)]; % etapa para deixar o tamanho do tmp2 igual ao vetor batimentos
        tmp2 = abs(tmp-a2); %valores da diferenca a partir do ponto máximo
        %t1
        [idy1 idx1] = min(tmp1); %indice do t1
        %t2
        [idy2 idx2] = min(tmp2); %indice t2
        tempos(1,i) = idx1;
        tempos(2,i) = idx2;
        tempos(3,i) = (idx1 + idx2)/2;
        alinhado_dn(:,i) = batimentos(tempos(1,i)-100:tempos(1,i)+100,i); 
    end
%Sinal médio final (DN)
medio_dn = 0;
for i = 1:400
   medio_dn = medio_dn + alinhado_dn(:,i);
end
medio_dn = medio_dn/400;


    % Integral Dupla
    D = zeros(1,400);
    modelo = s(T1-20:T1+60,1);
    for i = 1:400;
        lag = zeros(2,11);
        lag(1,:)=(-5:5);
        modeloid = cumtrapz((cumtrapz((modelo).^2)).^2);
        %modeloid = modelo;
        for l = -5:+5
            id = cumtrapz((cumtrapz((batimentos(T1-20+l:T1+60+l,i).^2))).^2);
            %id = batimentos(T1-20+l:T1+60+l,i);
            cc = corrcoef(id,modeloid);
            lag(2,l+6) = cc(1,2);
        end
        [c,ind] = max(lag(2,:));
        d = lag(1,ind);
        D(i) = d;
        alinhado_id(:,i) = batimentos(T1-100+d:T1+100+d,i);
    end
%Sinal médio final (ID)
medio_id = 0;
for i = 1:400
    medio_id = medio_id + alinhado_id(:,i);
end
medio_id = medio_id/400;

    
    %Integral Normalizada
    DD = zeros(1,400);
    mod = s(T1-20:T1+60,1);
    for i = 1:400;
        modtri = mod.^2;
        modni = cumtrapz(modtri)/max(cumtrapz(modtri));
        tri = (batimentos(T1-20:T1+60,i)).^2;
        ni = cumtrapz(tri)/max(cumtrapz(tri));
        d = trapz(modni - ni);
        d = round(d);
        DD(i) = d;
        alinhado_in(:,i) = batimentos(T1-100+d:T1+100+d,i);
    end
    
%Sinal médio final (IN)
medio_in = 0;
for i = 1:400
   medio_in = medio_in + alinhado_in(:,i);
end
medio_in = medio_in/400;

    
   % Integral Dupla com detector de nivel
   D0 =  zeros(1,400);
   D1 = zeros(1,400);
   modelo = s(T1-20:T1+60,1);
    for i = 1:400;
        %modelo = s(tempos(1,i)-20:tempos(1,i)+60,1); %linha adicionda
        lag = zeros(2,11);
        lag(1,:)=(-5:5);
        modeloid = cumtrapz((cumtrapz((modelo).^2)).^2);
        %modeloid = modelo;
        for l = -5:+5
            id = cumtrapz((cumtrapz((batimentos(tempos(1,i)-20+l:tempos(1,i)+60+l,i).^2))).^2);
            %id = batimentos(T1-20+l:T1+60+l,i);
            cc = corrcoef(id,modeloid);
            lag(2,l+6) = cc(1,2);
        end
        [c,ind] = max(lag(2,:));
        d = lag(1,ind);
        D0(i) = d;
        D1(i) = d + tempos(1,i);
    end
    
    %Integral Normalizada com detector de nivel
    DD0 = zeros(1,400); 
    DD1 = zeros(1,400);
    mod = s(T1-20:T1+60,1);
    for i = 1:400;
        %mod = s(tempos(1,i)-20:tempos(1,i)+60,1); %linha adicionda
        modtri = mod.^2;
        modni = cumtrapz(modtri)/max(cumtrapz(modtri));
        tri = (batimentos(tempos(1,i)-20:tempos(1,i)+60,i)).^2;
        ni = cumtrapz(tri)/max(cumtrapz(tri));
        d = trapz(modni - ni);
        %d = round(d);
        DD0(i) = d;
        DD1(i) = d + tempos(1,i);
    end
    dn = std(tempos(1,:));
    id = std(D0);
    in = std(D0);
    dnid = std(D1);
    dnin = std(DD1) ;
    
    DN(ii) = dn;
    ID(ii) = id;
    IN(ii) = in;
    DNID(ii) = dnid;
    DNIN(ii) = dnin;
    
    % COVARIANCIA
    %     CDNID = cov(tempos(1,:),D0);
    %     CDNNI = cov(tempos(1,:),DD0);
    %     CovDNID = CDNID(1,2)
    %     CovDNNI = CDNNI(1,2)
    
end
%DN
mDN = mean(DN)
stdDN = std(DN)

%ID
mID = mean(ID)
stdID = std(ID)

%IN
mIN = mean(IN)
stdIN = std(IN)

%DNID
mDNID = mean(DNID)
stdDNID = std(DNID)

%DNIN
mDNIN = mean(DNIN)
stdDNIN = std(DNIN)

%Detector de nivel
% figure
% hold on
% for i = 1: 400
%    plot((batimentos(tempos(1,i)-20:tempos(1,i)+60,i))) 
% 
% end
% mediadn = 0;
% for i = 1:400
%    mediadn = mediadn + batimentos(tempos(1,i)-20:tempos(1,i)+60,i);
% end
% mediadn = mediadn/400;
% plot(s(T1-20:T1+60,1),'r')
% plot(mediadn,'y')
% com o aumento do ruido a moda de tempos(1,:) vai divergindo de T1;

%integral dupla
% figure
% hold on
% for i = 1: 400
%    plot((batimentos(T1+D(i)-20:T1+D(i)+60,i))) 
%    
% end
% mediaid = 0;
% for i = 1:400
%    mediaid = mediaid + batimentos(T1+D(i)-20:T1+D(i)+60,i);
% end
% mediaid = mediaid/400;
% plot(s(T1-20:T1+60,1),'r')
% plot(mediaid,'y')

% DN = std(tempos(1,:))
% %ID = std(D) %precisa ser o desvio padrão do D0
% ID = std(D0)
% %NI = std(DD) %precisa ser o desvio padrão do DD0
% IN = std(D0)
% DNID = std(D1)
% DNIN = std(DD1) 
% 
% CDNID = cov(tempos(1,:),D0);
% CDNNI = cov(tempos(1,:),DD0);
% CovDNID = CDNID(1,2)
% CovDNNI = CDNNI(1,2)

% 
% % Espectro de frequencia de um batimento
% sbat = batimentosf(T2-70:T2+70,1);
% n1 = size(sbat);
% fa1 = 1000;
% L1 = length(sbat');
% La = round(L1/2);
% winvec1 = hann(L1); 
% F = fa1/L1*(0:La);
% x = sbat;
% X = fft(x.*winvec1);
% Px1 = X.*conj(X)/L1;
% dBx = pow2db(Px1);
% figure
% plot(F(1:round(La)),dBx(1:round(La)))
% title('Espectro de um batimento', 'fontsize', 18);
% xlabel('Frenquência (Hz)','fontsize',18)
% ylabel('Potência', 'fontsize',18)
% set(gca,'fontsize',16);
% 
% % Espectro de frequencia dos batimentos promediados
% [~,T3] = max(medio_id);
% n1 = size(medio_id(T3-70:T3+70));
% fa1 = 1000;
% L1 = length(medio_dn');
% La = round(L1/2);
% winvec1 = hann(L1); 
% F = fa1/L1*(0:La);
% y1 = medio_dn;
% Y1 = fft(y1.*winvec1);
% Py1 = Y1.*conj(Y1)/L1;
% dBy = pow2db(Py1);
% 
% figure
% plot(F(1:round(La)),dBy(1:round(La)))
% title('Espectro dos batimentos promediados', 'fontsize', 18);
% xlabel('Frenquência (Hz)','fontsize',18)
% ylabel('Potência', 'fontsize',18)
% set(gca,'fontsize',16);
% 
% P = (Py1(1:round(La))./(Px1(1:round(La))));
% dBp = pow2db(P);
% figure
% plot(F(1:round(La)),dBp)
% title('Função de transferencia do ECGAR', 'fontsize', 18);
% xlabel('Frenquência (Hz)','fontsize',18)
% ylabel('Potência [dB]', 'fontsize',18)
% set(gca,'fontsize',16);
% grid on
% Plotar em semilog!!!!!!!!!!!!!!!!!!!!!!

% % Estimando a função de transferência
% Y = medio_dn(76:166);
% X = batimentosf(T1-54:T1+36,4);
% fa = 1000;
% figure
% tfestimate(X,Y,8,[],[],fa)