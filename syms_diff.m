% ######################################################################## 
%                        LEZIONE SISTEMI UNDERDAMPED
% Analisi di sistemi armonici soggetti ad impulso e con smorzamento, e
% risposta in frequenza di tali sistemi
% ########################################################################
% 
% Tale file va eseguito una tantum per ottenere le equazioni del moto, che
% possiamo poi ricopiare dalla command window e inserire nel nostro script
% 
% ########################################################################
%
% Filippo Badalamenti, 22/11/19

clear all
clc

%% Ottengo le equazioni

% variabii rese simboliche per la derivata
syms A zeta wn wd phi t

% espressione della legge oraria nel tempo di un sistema del 2° ordine
z_t = A*exp(-zeta*wn*t).*sin(wd*t+phi);

% derivo simbolicamente
v_t = diff(z_t, t) %derivo rispetto al simbolo t

% derivo simbolicamente
a_t = diff(v_t, t)