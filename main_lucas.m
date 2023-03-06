%%  S5 - APP4 - PROBLEMATIQUE - MAIN.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Lucas Corrales
%   CIP:        CORL0701

%   Date:       2-MARS-2023
%   Modifications (Date - initiales - détails):


% JAMAIS UN BODE AVEC UNE FTBF!!!!!!!

%% MAIN CODE
clc
close all
clear all

% Initialisation
constantes % call le fichier des constantes



%% DEBUT DE La problematique

%ESSAIE #1
%syms s
%s=tf('s');

%I = eye(4,4);

%inter = (s.*I)-A;
%inv_inter = inv(inter);

%TF = C*(inv_inter*B);




% ESSAIE #2
% Compute the transfer function
%[num,den] = ss2tf(A,B,C,D);

% Print the transfer function
%tf = tf(num,den)
%disp(tf)









%% QUESTION a)

disp(['------QUESTION a)------']);

% Calcul des valeurs propres du systeme
val_propres_A = eig(A)

% calcul des carac temporelles
wn = abs(val_propres_A);
zeta = -real(val_propres_A)./wn;
wa = wn.*sqrt(1-zeta.^2);
phi = acos(zeta);
Mp = 100*exp(-pi./tan(phi));
ts = 4./(zeta.*wn);
tp = pi./wa;

% affichage des carac temporelles
disp(["Affichage des carac temporelles:"]);
disp(['wn = ', num2str(wn(end)), ' rad/s']);
disp(['zeta = ', num2str(zeta(end)), ' unites']);
disp(['wa = ', num2str(wa(end)), ' rad/s']);
disp(['phi = ', num2str(phi(end)), ' radian']);
disp(['Mp = ', num2str(Mp(end)), ' %']);
disp(['ts = ', num2str(ts(end)), ' s']);
disp(['tp = ', num2str(tp(end)), ' s']);






%% QUESTION b) - v/a_prop
[num,den] = ss2tf(A,B,C,D,2);   % le 2 signifie quon veut le 2e element de U soit aprop
v_aprop = tf(num(1, :), den)         % v/aprop
alpha_aprop = tf(num(2, :), den);         % alpha/aprop
teta_aprop = tf(num(3, :), den);         % teta / aprop
q_aprop = tf(num(4, :), den);         % q/aprop
gamma_aprop = tf(num(5, :), den);         % gamma/aprop


numerateur = num(1,:);
zeros = roots(num(1,:))
poles = roots(den)


p1 = poles(1);
p2 = poles(2);
p3 = poles(3);
p4 = poles(4);


z1 = zeros(1);
z2 = zeros(2);
z3 = zeros(3);


figure('Name', 'Rlocus de V/a_prop')
rlocus(v_aprop)

figure('Name', 'Step de V/a_prop')
step(v_aprop)


[num, den] = zp2tf(zeros, poles, 1);
TF_a_partir_zeros_poles_v_aprop = tf(num,den)
 

% section lucas
% identification de la fonction de transfert à partir des pôles et zéros
z = zeros;
p = val_propres_A.';
%p = poles;
k = 1;
sys_zpk = zpk(z, p, k);
display(sys_zpk)






%% QUESTION b) - phase non minimale
disp(['------QUESTION b)------']);
[num,den] = ss2tf(A,B,C,D,1);   % le 1 signifie quon veut le 1er element de U soit delta_c
v_delta_c = tf(num(1, :), den)         % v/delta_c

figure('Name', 'V sur Delta_C')
rlocus(v_delta_c)

figure('Name', 'V sur Delta_C')
bode(v_delta_c), grid










%% FEEDBACK d)
disp(['------QUESTION d)------']);

% EXEMPLE FEEDBACK
%num = [1 2];
%den = [1 3 2];
%G = tf(num, den)
%Kv = 0.5;
%H = tf(Kv, 1);
%T = feedback(G*H, 1);
%step(T)


G = v_aprop
% Kv = [0     0.1     0.5     1     25  50  75  100       200     300     500     1000];
% Kv = [0     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9     50000];
Kv = 1.38

figure('Name', 'Kv')
for i = 1:length(Kv)
    H = tf(Kv(i), 1);
    T = feedback(G*H, 1);
%     step(T)
    rlocus(T)
    legend
    hold on
end
[num_T,den_T] = tfdata(T)
[A1,B1,C1,D1] = tf2ss(num_T{1},den_T{1})





%% question e)
disp(['------QUESTION e)------']);

% choix de Kv
% Calcul des valeurs propres du systeme
intermediaire = v_aprop.*Kv
val_propres_v_aprop = eig(intermediaire)


% calcul des carac temporelles
wn = abs(intermediaire);
zeta = -real(intermediaire)./wn;
wa = wn.*sqrt(1-zeta.^2);
phi = acos(zeta);
Mp = 100*exp(-pi./tan(phi));
ts = 4./(zeta.*wn);
tp = pi./wa;

% affichage des carac temporelles
disp(["Affichage des carac temporelles:"]);
disp(['wn = ', num2str(wn(end)), ' rad/s']);
disp(['zeta = ', num2str(zeta(end)), ' unites']);
disp(['wa = ', num2str(wa(end)), ' rad/s']);
disp(['phi = ', num2str(phi(end)), ' radian']);
disp(['Mp = ', num2str(Mp(end)), ' %']);
disp(['ts = ', num2str(ts(end)), ' s']);
disp(['tp = ', num2str(tp(end)), ' s']);




% calcul du nouveau modele n1_d1
%num = [wn.^2]
%den = [1  ,     2.*zeta.*wn,       wn.^2]

%n1_d1 = tf(num, den)




%% question f)
disp(['------QUESTION f)------']);
% 

% Kv = [0     0.1     0.5     1     25  50  75  100       200     300     500     1000];
% Kv = [0     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9     50000];
Kv = 1.386

figure('Name', 'Bode')
for i = 1:length(Kv)  
    bode(Kv(i)*G), grid
    margin(Kv(i)*G);
    hold on
end

fprintf('Phase margin: %0.2f degrees\n', Pm);
fprintf('Gain margin: %0.2f dB\n', Gm);
fprintf('Omega G: %0.2f rad/s\n', Wog);
fprintf('Omega P: %0.2f rad/s\n', wop);

%% question f)
%trouver les pôles dominant (les pôle qui affect le plus la fonction)
[R,P,K]=residue(G.numerator{1},G.denominator{1});

%trouver le poid des pôles
Cdom=abs(R)./abs(real(P))
 
%le poid pole avec la plus grande valeur=plus proche de l'axe imaginaire
%on veut toujours deux pole... systeme d'ordre 2

%reduction du systeme
[num,den]=residue(R(3:4),P(3:4),K);
 
TFR=tf(num,den);

g0 = dcgain(G);
g1 = dcgain(TFR);

G_simpl = (g0/g1)*TFR

step(G)
hold on
step(G_simpl)
hold on
step(v_aprop)









