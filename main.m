%%  S5 - APP4 - PROBLEMATIQUE - MAIN.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Lucas Corrales
%   CIP:        CORL0701

%   Date:       2-MARS-2023
%   Modifications (Date - initiales - détails):


%% MAIN CODE
clc
close all
clear all

% Initialisation
constantes % call le fichier des constantes







%% ANNEXE A

% instanciation des matrices 
A = [   -0.018223   -0.088571  -9.78   0;...
        -0.003038    -1.2563     0     1;...
            0            0       0     1;...
         0.0617       -28.075    0   -4.5937];


B = [     0        1.1962;...
          0       -0.00120;...
          0           0;...
          7.84      -4.05];


C = [     1        0        0        0;...
          0       57.296    0        0;...
          0        0      57.296     0;...
          0        0        0    57.296;...
          0      -57.296    57.296   0];


D = [   0   0;...
        0   0;...
        0   0;...
        0   0;...
        0   0];


% Variables d'entree
U = [delta_c   a_prop]' ;     %en degres et en fraction de la poussee maximale


% Variables d'etat
X = [v   alpha   teta   q]';     % en m/s et en radians


% Variables de sortie
Y = [v  alpha   teta    q   gamma]';     %en m/s et en degres











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









% ESSAIE #3
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
disp(["Affichage des carac temporelles:"])
disp(['wn = ', num2str(wn(end)), ' rad/s']);
disp(['zeta = ', num2str(zeta(end)), ' unites']);
disp(['wa = ', num2str(wa(end)), ' rad/s']);
disp(['phi = ', num2str(phi(end)), ' radian']);
disp(['Mp = ', num2str(Mp(end)), ' %']);
disp(['ts = ', num2str(ts(end)), ' s']);
disp(['tp = ', num2str(tp(end)), ' s']);






%% FONCTIONS DE TRANSFERT
% Compute the transfer function
[num,den] = ss2tf(A,B,C,D,2);   % le 2 signifie quon veut le 2e element de U soit aprop
v_aprop = tf(num(1, :), den)         % v/aprop
alpha_aprop = tf(num(2, :), den)         % alpha/aprop
teta_aprop = tf(num(3, :), den)         % teta / aprop
q_aprop = tf(num(4, :), den)         % q/aprop
gamma_aprop = tf(num(5, :), den)         % gamma/aprop


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


figure(1)
rlocus(v_aprop)

figure(2)
step(v_aprop)


[num, den] = zp2tf(zeros, poles, 1);
TF_non_min = tf(num,den)


%% section lucas
% identification de la fonction de transfert à partir des pôles et zéros
z =zeros;
p = val_propres_A.';
%p = poles;
k = 1;
sys_zpk = zpk(z, p, k);
display(sys_zpk)


%% FEEDBACK

% EXEMPLE FEEDBACK
%num = [1 2];
%den = [1 3 2];
%G = tf(num, den)
%Kv = 0.5;
%H = tf(Kv, 1);
%T = feedback(G*H, 1);
%step(T)


G = TF_non_min
Kv = [0     0.1     0.5     1     25  50  75  100       200     300     500     1000];
Kv = [0     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9     50000];


for i = 1:length(Kv)
    H = tf(Kv(i), 1);
    T = feedback(G*H, 1);
    step(T)
    legend
    hold on
end








