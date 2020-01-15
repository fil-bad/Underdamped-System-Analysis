% ######################################################################## 
%                        LEZIONE SISTEMI UNDERDAMPED
% Analisi di sistemi armonici soggetti ad impulso e con smorzamento, e
% risposta in frequenza di tali sistemi
% ########################################################################
% Scaletta:
% 1) Caricamento dati ottenuti da file esterno '.txt'
% 2) Smoothing dati acquisiti
% 3) Analisi frequenziale tramite FFT
% 4) Analisi di una singola componente frequenziale, e ricavare la sua
%    legge oraria tramite Logarithmic Decrement e fit con esponenziale
% 5) Ottenimento del sistema di Laplace, per la simulazione con impulso 
% 
% ########################################################################
% 
% Filippo Badalamenti, 22/11/19

clear all
clc
close all
%% Caricamento misure da file esterno

% filename = 'Misure1brevesyncImp.txt';
filename = 'Test1vecchio.txt';

delimiterIn = '\t'; % separatore dei dati raccolti, per gli assi [x y z]
Measures = importdata(filename,delimiterIn)*9.81; %% portiamo in m/s^2 
% matrice Nx3, con N misure effettuate

%% Inizializzazione variabili usate globalmente nel programma
% NB sto misurando le accelerazioni!!! non dovrà coincidere
% il risultato della curva con la sua equazione del moto; basta
% però ricordare che a(t) ~= -(wd^2)*x(t)

x = Measures(:,1);
y = Measures(:,2);
z = Measures(:,3);

len = length(z);

Ts = 0.02; % tempo di campionamento Arduino
Fs = 1/Ts;
t = [0:Ts:(len-1)*Ts]'; % vettore dei tempi


figure(1)
plot(t, x, 'r', t, y, 'g', t, z, 'b')
legend('x','y','z')

%% Curve smoothing

t_dense=[0:Ts/10:(len-1)*(Ts)]';
% vado ad interpolare i dati raccolti tramite spline (curva morbida)
z_dense = interp1(t,z,t_dense,'spline');

figure(2)
plot(t,z,'o',t_dense,z_dense,':.');
legend('z', 'z_d_e_n_s_e')

% t = t_dense;
Ts_d = Ts/10;
Fs_d = 1/Ts_d;
len_d = length(z_dense);
%z = z_dense;

%% Analisi frequenziale Fourier

% Da qui in poi analizziamo solo z, la procedura è ripetibile 
% per qualunque asse

Z_f = fft(z_dense); % analisi frequenze di z_smoothed
P2 = abs(Z_f/len_d);
P1 = P2(1:len_d/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs_d*(0:(len_d/2))/len_d;

zfil = bandpass(z,[2 (Fs/2)-2],Fs); 
% l'analisi sulle accelerazioni funziona bene da frequenza ca 10^2
% l'accelerazione è un passa alto -> tolgo le frequenze basse
Z_fFil = fft(zfil);
P2_z = abs(Z_fFil/len);
P1_z = P2_z(1:len/2+1);
P1_z(2:end-1) = 2*P1_z(2:end-1);
f2 = Fs*(0:(len/2))/len;

figure(3)
plot(f,P1,'b')%, f2,P1_z,'r.')
legend('P1_z','P1_z_f_i_l')
title('Single-Sided Amplitude Spectrum of Z_d_e_n_s_e(t) & Z(t)')
xlabel('f(Hz)')
ylabel('|P1(f)|,|P1_z(f)|')
xlim([0 30]) % per centrare la vista su un intervallo delle frequenze

%% Analizzo d'ora in poi il segnale smoothed

z_old = z; % salvo i vecchi valori, mi serviranno per il modello ARMA
len_old = len;
ts_old = Ts;
t_old = t;

t = t_dense;
Ts = Ts_d;
Fs = Fs_d;
len = len_d;
z = z_dense;

%% Ottengo i maggiori contributi frequenziali
% prendo i picchi dal grafico di Fourier, considerando quelli che possano
% contribuire fino al 20% del segnale

[ampl, locs] = findpeaks(P1,'MinPeakHeight',max(P1)/5,'MinPeakDistance',1);

freqs = zeros(length(ampl),1);

for i=1:1:length(ampl)

    freqs(i) = f(locs(i));

end

%% Segnale ricostruito prendendo solo alcune componenti

z_tot = zeros(len,1);

for k=1:1:length(freqs)
    
    z_tot = z_tot + bandpass(z,[freqs(k)-1 freqs(k)+1],Fs);
end

figure(4)
plot(t,z_tot,'b', t,z,'r')
legend('z_t_o_t','z_a_c_q_u_i_s_i_t_o')

%% Visualizzo solo la frequenza di interesse

% seleziono solo un piccolo intervallo intorno al picco d'interesse
z_fil = bandpass(z,[freqs(1)-1 freqs(1)+1],Fs); 
z_Ffil = fft(z_fil);
P2_zf = abs(z_Ffil/len);
P1_zf = P2_zf(1:len/2+1);
P1_zf(2:end-1) = 2*P1_zf(2:end-1);
f_zf = Fs*(0:(len/2))/len;

figure(5)
plot(f,P1,'r.', f_zf,P1_zf,'b')
legend('P1_z','P1_zf')
title('Single-Sided Amplitude Spectrum of Z(t)')
xlabel('f(Hz)')
ylabel('|P1_z(f)|')
xlim([0 30])

%% Calcolo zeta e wn tramite Peak-Picking e Damping Ratio

[amp, peaks] = findpeaks(z_fil); 

figure(4)
plot(t,z_fil,'b',peaks*Ts,amp,'r*');

% wd_bad = mean(diff(peaks)) % far vedere su command window

wd = freqs(1)*2*pi; %% meglio questa per l'analisi di wd;
% fourier è più preciso rispetto alla media dei picchi


fitted = fit(peaks,amp,'exp1')

A = fitted.a; % trovo l'ampiezza A

% prendo due picchi a caso per il damping ratio;
% trovando l'esponenziale in realtà, potevo fare 
% zeta*wn = fitted.b , e risolvendo le equazioni

n1 = 5;
n2 = 15;

sigma = (1/(n2-n1))*log(z_fil(peaks(n1))/z_fil(peaks(n2)));

zeta = 1/sqrt(1+((2*pi)/sigma)^2 );

wn = wd/sqrt(1-zeta^2);

phi = atan2(sqrt(1-zeta^2), zeta);

% questo perchè la legge oraria è:
% z1_t = ( A*exp(-zeta*wn*t).*sin(wd*t+phi) );


%% Ricavo le equazioni del moto

z1_t = ( A*exp(-zeta*wn*t).*sin(wd*t+phi) );

% equazione a1_t ottenuta tramite syms_diff.m
a1_t = A*wn^2*zeta^2*exp(-t*wn*zeta).*sin(phi + t*wd) - A*wd^2*exp(-t*wn*zeta).*sin(phi + t*wd) - 2*A*wd*wn*zeta*exp(-t*wn*zeta).*cos(phi + t*wd);
% N.B.: noi abbiamo preso solamente le accelerazioni! quindi la curva che
% abbiamo trovato non è altro che a1_t; poichè lo smorzamento è limitato,
% possiamo approssimare alla relazione del moto armonico semplice,
% ovvero a(t) ~= -(wd^2)*x(t)

figure(10)
plot(t, z_fil,'b', t,z1_t,'r')% -(1/wn^2)*a1_t,'r')
legend('z_f_i_l(t)', 'z1(t)')


%% Simulo il sistema differenziale

z_fil_old = bandpass(z_old,[freqs(1)-1 freqs(1)+1],1/ts_old); 

% ricostruisco il segnale con i tool di Matlab di identificazione
tf_z1 = ssest(z_fil_old, 2, 'Ts', ts_old)

tf_c = d2c(tf_z1) % ottengo il sistema continuo

[wn_v2,zeta_v2,p2] = damp(tf_c.A); % ottengo i dati dello stato

% riottengo i parametri del sistema, per vederlo nel tempo
wn2 = wn_v2(1)
zeta2 = zeta_v2(1)*2 % necessario, per lo smorzamento
wd2 = wn2*sqrt(1-zeta2^2);
phi2 = atan2(sqrt(1-zeta2^2), zeta2)

z1_t_ssest = ( A*exp(-zeta2*wn2*t_old).*sin(wd2*t_old+phi2) );
figure(15)
plot(t_old, z_fil_old,'b', t_old,z1_t_ssest,'r')% -(1/wn^2)*a1_t,'r')
legend('z_f_i_l(t)', 'z1_s_s_e_s_t(t)')


%% Riscrivo il sistema in Laplace

% ottengo il polinomio caratteristico del sistema, in modo da inserirlo
% nella funzione di trasferimento

pol_car = charpoly(tf_c.A);
pol_car(2) = 2*pol_car(2);

% sistema ottenuto con i tool di automatici
sys_c = tf(wn2^2, pol_car)
% sistema ottenuto con i metodi del prof, più preciso
sys_rec = tf(wn^2, [1 2*zeta*wn wn^2])  

% Ho riscritto tali funzioni perchè è un sistema del secondo ordine
% smorzato, e che ha la sua forma standard pari :
%
% sys = wn^2/(s^2 + s*(2*zeta*wn) + wn^2)



% N.B. : Abbiamo ottenuto così un polinomio di Laplace, che permette di
% lavorare con numerosi strumenti studiati in Automazione. La cosa
% positiva? possiamo risalire alla risposta in velocità e in posizione
% moltiplicando anzichè integrando (vantaggio di Laplace)

%% Simulo la risposta all'impulso dell'accelerazione

figure(20)
impulse(sys_c, 'r')
hold on
impulse(sys_rec, 'b')
legend('sys_c', 'sys_rec')


