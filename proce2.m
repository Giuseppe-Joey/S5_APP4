%Q13
%a)
G = tf([17.8885], [1 2 0])
figure()
bode(G),grid
hold on
bode(G*Ga),grid
figure()
subplot(211)
rlocus(G)
subplot(212)
rlocus(G*Ga)

%b)
Ga =1.5697 * tf([1 2.5483], [1 6.2787])
figure()
bode(Ga),grid

%e)
Pd = 


[Gm,Pm,Wp,Wg] = margin(MAG,PHASE,W)