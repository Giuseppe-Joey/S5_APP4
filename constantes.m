%%  S5 - APP4 - PROBLEMATIQUE - CONSTANTES.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Lucas Corrales
%   CIP:        CORL0701

%   Date:       2-MARS-2023
%   Modifications (Date - initiales - détails):

% Description: fichier contenant toutes les constantes necessaires 


%% ANNEXE A

delta_c = 1;    % angle du gouvernail de profondeur
a_prop = 0.5;     %fraction de la poussee maximale des moteurs (0 < a_prop < 1)
v = 23.4;       %vitesse de lavion
alpha = 30;     %angle dattaque(entre laxe longitudinal de lavion et la vitesse)


teta = 15;      %angle de tangage(entre laxe longitudinal de lavion et le plan horizontal)

q = 45;         %vitesse angulaire en tangage de lavion (q=d(teta)/dt)

gamma = teta - alpha;     %angle de vol de lavion (gamma = teta - alpha) et (gamma > 0:lavion monte, gamma < 0:lavion descend)








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




