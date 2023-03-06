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
propres = eig(A)

% calcul des carac temporelles
wn = abs(propres);
zeta = -real(propres)./wn;
wa = wn.*sqrt(1-zeta.^2);
phi = acos(zeta);
Mp = 100*exp(-pi./tan(phi));
ts = 4./(zeta.*wn);
tp = pi./wa;

% affichage des carac temporelles
disp(["Valeurs propres du systeme:"])
disp(['tp = ', num2str(tp(end)), ' (secs)'])
disp(['ts = ', num2str(ts(end)), ' (secs)'])



