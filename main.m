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

disp(['------Initialisation------']);

% Initialisation
constantes % call le fichier des constantes

disp(['']);
disp(['']);





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


disp(['']);
disp(['']);






%% QUESTION b) - v/a_prop
disp(['------QUESTION b) - v/aprop------']);
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


figure('Name', 'Rlocus de V/a_prop');
rlocus(v_aprop);

figure('Name', 'Step de V/a_prop');
step(v_aprop);


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

disp(['']);
disp(['']);




%% QUESTION b) - phase non minimale
disp(['------QUESTION b)------']);
[num,den] = ss2tf(A,B,C,D,1);   % le 1 signifie quon veut le 1er element de U soit delta_c
v_delta_c = tf(num(1, :), den)         % v/delta_c

figure('Name', 'V sur Delta_C');
rlocus(v_delta_c);

figure('Name', 'V sur Delta_C');
bode(v_delta_c), grid;

disp(['']);
disp(['']);








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


G = TF_a_partir_zeros_poles_v_aprop
Kv = [0     0.1     0.5     1     25  50  75  100       200     300     500     1000];
Kv = [0     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9     50000];

figure('Name', 'Kv');
for i = 1:length(Kv)
    H = tf(Kv(i), 1);
    T = feedback(G*H, 1);
    step(T);
    legend
    hold on
end


disp(['']);
disp(['']);








%% question e)
disp(['------QUESTION e)------']);
% Kv = [0     0.1     0.5     1     25  50  75  100       200     300     500     1000];
% Kv = [0     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9     50000];

G = v_aprop
Kv = 1.0262;

    figure('Name', 'Kv');
    H = tf(Kv,1);
    T = feedback(G*H, 1);
    rlocus(T);
    legend

[num_T,den_T] = tfdata(T)
[A1,B1,C1,D1] = tf2ss(num_T{1},den_T{1})


disp(['']);
disp(['']);





%% question f)
disp(['------QUESTION f)------']);
[Gm, Pm, Wog, wop] = margin(G);

fprintf('Phase margin: %0.2f degrees\n', Pm)
fprintf('Gain margin: %0.2f dB\n', Gm)
fprintf('Omega G: %0.2f rad/s\n', Wog)
fprintf('Omega P: %0.2f rad/s\n', wop)

figure('Name', 'Bode');
bode(G), grid;


disp(['']);
disp(['']);








%% question g)
disp(['------QUESTION g)------']);

%trouver les pôles dominant (les pôle qui affect le plus la fonction)
[R,P,K]=residue(G.numerator{1},G.denominator{1});

%trouver le poid des pôles
Cdom=abs(R)./abs(real(P))
 
%le poid pole avec la plus grande valeur=plus proche de l'axe imaginaire
%on veut toujours deux pole... systeme d'ordre 2
%le pole 3 et 4 sont ceux qui ont le plus de poid

%reduction du systeme
[num,den]=residue(R(3:4),P(3:4),K);
 
TFR=tf(num,den);

g0 = dcgain(G);
g1 = dcgain(TFR);

G_simpl = (g0/g1)*TFR

[num_G_simpl,den_G_simpl] = tfdata(G_simpl)

Wn_G_simpl = sqrt(den_G_simpl{1}(3))
Zeta = (den_G_simpl{1}(2))/(2*Wn_G_simpl)

Kv1 = ((Wn_G_simpl^2)-0.047)/1.465
Kv2 = ((den_G_simpl{1}(2))-0.028)/2.122
 
figure('Name', 'Comparaison avec nouvelle tf');
step(G);
hold on
step(G_simpl);
legend



disp(['']);
disp(['']);









%% question h)
disp(['------QUESTION h)------']);




disp(['']);
disp(['']);






%% question i)
disp(['------QUESTION i)------']);
%trouver les pôles dominant (les pôle qui affect le plus la fonction)
[R,P,K]=residue(v_aprop.numerator{1},v_aprop.denominator{1});

%trouver le poid des pôles
Cdom=abs(R)./abs(real(P));
 
%le poid pole avec la plus grande valeur=plus proche de l'axe imaginaire
%on veut toujours deux pole... systeme d'ordre 2

%reduction du systeme
[num,den]=residue(R(3:4),P(3:4),K);
 
FT_originale = v_aprop
FT_reduite = tf(num,den)


figure;
rlocus(FT_reduite), grid
hold on
rlocus(FT_originale), grid
legend


disp(['']);
disp(['']);





%% question j)
disp(['------QUESTION j)------']);

% Calcul des valeurs propres du systeme
%val_propres_A1 = eig(A1)

% [num,den] = ss2tf(A1,B1,C1,D1,1);   % le 1 signifie quon veut le 1e element de U soit delta
% gamma_delta = tf(num(5, :), den)         % gamma/delta
% 
% figure('Name', 'Rlocus de gamma/delta');
% rlocus(gamma_delta);
% 
% figure('Name', 'Step de gamma/delta');
% step(gamma_delta);
% 
% 
% % calcul des carac temporelles
% wn = abs(val_propres_A1);
% zeta = -real(val_propres_A1)./wn;
% wa = wn.*sqrt(1-zeta.^2);
% phi = acos(zeta);
% Mp = 100*exp(-pi./tan(phi));
% ts = 4./(zeta.*wn);
% tp = pi./wa;
% 
% 
% % affichage des carac temporelles
% disp(["Affichage des carac temporelles:"]);
% disp(['wn = ', num2str(wn(end)), ' rad/s']);
% disp(['zeta = ', num2str(zeta(end)), ' unites']);
% disp(['wa = ', num2str(wa(end)), ' rad/s']);
% disp(['phi = ', num2str(phi(end)), ' radian']);
% disp(['Mp = ', num2str(Mp(end)), ' %']);
% disp(['ts = ', num2str(ts(end)), ' s']);
% disp(['tp = ', num2str(tp(end)), ' s']);


disp(['']);
disp(['']);








%% question k)
disp(['------QUESTION k)------']);





disp(['']);
disp(['']);




%% question l)
disp(['------QUESTION l)------']);








disp(['']);
disp(['']);



%% question m)
disp(['------QUESTION m)------']);

