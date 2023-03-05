clc
clear 
close all

%probleme 6
%6(a)
%on veut en FTBO
wn=1;
zeta=0.1;
kp=1;
nr=[wn^2];
dr=[1 2*wn*zeta wn^2]
FT=tf(nr,dr)

td=[0.0 0.2 0.5 1.0 2.5 5.0];

figure(1)
figure(2)
for i=1:6
num=[kp*td(i) kp];
den=1;
FT2=tf(num,den);
FTBO=FT*FT2;
figure(1)
%subplot(2,1,1)

rlocus(FTBO)
hold on
p=rlocus(FTBO,1);%calcul les poles de la ftbf a k=1 
plot (real(p), imag(p),'p')

figure(2)
%subplot(2,1,2)
t=[0:0.01:50]';
u=ones(size(t));
FTBF=feedback(kp*FTBO,1)
lsim(FTBF,u,t)


%    FTBF=feedback(kp*FT,1);
%    DEN=[kp*td(i) kp];
%    TF=tf(NUM,DEN);
%    lsim(TF,u,t);
    hold on
end

%%
clc
clear 
close all

%probleme 6
%6(b)
%on veut en FTBO
wn=1;
zeta=0.1;
kp=1;
nr=[wn^2];
dr=[1 2*wn*zeta wn^2];
FT=tf(nr,dr)

tau=[0.0 0.01 0.05 0.1 0.2];
figure(1)
figure(2)
for i=1:5
num=[1];
den=[tau(i) 1];
FT2=tf(num,den);
FTBO=FT*FT2;
figure(1)
%subplot(2,1,1)

rlocus(FTBO)
hold on
p=rlocus(FTBO,1);%calcul les poles de la ftbf a k=1 
plot (real(p), imag(p),'p')

figure(2)
%subplot(2,1,2)
t=[0:0.01:50]';
u=ones(size(t));
FTBF=feedback(kp*FTBO,1)
lsim(FTBF,u,t)


%    FTBF=feedback(kp*FT,1);
%    DEN=[kp*td(i) kp];
%    TF=tf(NUM,DEN);
%    lsim(TF,u,t);
    hold on
end
%%
%Q7

%probleme 7
%on veut en FTBO
clc
clear 
close all
wn=1;
zeta=0.1;
kp=1;
nr=[wn^2];
dr=[1 2*wn*zeta wn^2];
FT=tf(nr,dr)

tr=[0.05 0.1 0.2 0.25];
figure(1)
figure(2)
for i=1:4
    
[num,den]=pade(tr(i),4)
FT2=tf(num,den);
FTBO=FT*FT2;
figure(1)
%subplot(2,1,1)

rlocus(FTBO)
hold on
p=rlocus(FTBO,1);%calcul les poles de la ftbf a k=1 
plot (real(p), imag(p),'p')

figure(2)
%subplot(2,1,2)
t=[0:0.01:50]';
u=ones(size(t));
FTBF=feedback(kp*FTBO,1);
lsim(FTBF,u,t)
hold on
end

%%
%Q8
%8(a)
gc1=1;
gc2=4.68.*tf([1 2.9],[1 5.4])
gc3=4.68.*tf([1 5.4],[5.4])
FTp=tf(4,[1 2 0])
FT=gc1*FTp
[num,den]=tfdata(FT)
wn=sqrt(num{1}(3))
zeta=(den{1}(2))/(2*wn)
ts=4/(wn*zeta)
%%
%8b
%c<est un ordre 3
gc2=4.68.*tf([1 2.9],[1 5.4])
FTp=tf(4,[1 2 0])
FTBO=gc2*FTp
%[num,den]=tfdata(FT)
%wn=sqrt(num{1}(3))
%zeta=(den{1}(2))/(2*wn)
%ts=4/(wn*zeta)
%voire lecture chapitre 5 pytagore 
figure(1)
p=rlocus(FTBO,1)
rlocus(FTBO)
hold on
plot (real(p), imag(p),'p')
figure(2)
t=[0:0.01:5]';
u=ones(size(t));
FTBF=feedback(FTBO,1);
lsim(FTBF,u,t)

wn=abs(p)
zeta=-real(p)./wn
ts=4./(zeta.*wn)
%%
%8c
%pour gc1
gc1=1;
FTp1=tf(4,[1 2 0])
FTBO1=gc1*FTp

figure(1)
p1=rlocus(FTBO1,1)
rlocus(FTBO1)
hold on
plot (real(p1), imag(p1),'p')
figure(2)
t=[0:0.01:5]';
u=ones(size(t));
FTBF1=feedback(FTBO1,1);
lsim(FTBF1,u,t)
hold on
%pour gc2
gc2=4.68.*tf([1 2.9],[1 5.4])
FTp2=tf(4,[1 2 0])
FTBO2=gc2*FTp1

figure(1)
p2=rlocus(FTBO2,1)
rlocus(FTBO2)
hold on
plot (real(p2), imag(p2),'p')
figure(2)
FTBF2=feedback(FTBO2,1);
lsim(FTBF2,u,t)
hold on

%pour gc3
gc3=4.68.*tf([1 2.9],[5.4])
FTp=tf(4,[1 2 0])
FTBO=gc3*FTp

figure(1)
p=rlocus(FTBO,1)
rlocus(FTBO)
hold on
plot (real(p), imag(p),'p')
figure(2)

FTBF=feedback(FTBO,1);
lsim(FTBF,u,t)
%%
clc 
clear 
close all
%Q9
%9a analytiquement
%9bpt1
tr=2
[num,den]=pade(tr,2)
FTp1=tf(num,den)
%9bpt2 anlaytiquement
FTp2=tf(1,[1 3])
FTBO=FTp1*FTp2
%9c
[NUM,DEN]=pade(tr,3)
FTP=tf(NUM,DEN)
FTBO1=FTP*FTp2

figure(1)
%ordre 2
p=rlocus(FTBO,1)
rlocus(FTBO)
legend('ordre 1')
hold on
plot (real(p), imag(p),'p')
%ordre 3
p1=rlocus(FTBO1,1)
rlocus(FTBO1)
hold on
plot (real(p1), imag(p1),'p')
legend('ordre 2','pôles','ordre 3')
hold off

%ordre 2
t=[0:0.01:10]';
u=ones(size(t));

FTBF=feedback(FTBO,1);
%ordre 3
FTBF1=feedback(FTBO1,1);
figure(2)
lsim(FTBF,u,t)
hold on
lsim(FTBF1,u,t)
legend('ordre 2','ordre 3')

%%
%Q10a
clc
clear
close all
FTBO=tf(1,[1 8 19 12 0])
FTBF=feedback(FTBO,1);
p=rlocus(FTBO,1);
t=[0:0.01:10]';
u=ones(size(t));
figure(1)
rlocus(FTBO)
hold on
plot(real(p), imag(p),'p')
title('lieux des racine')
figure(2)
lsim(FTBF,u,t)
title('reponse à léchelon')
%%
%Q10b
clc
clear
close all
FTBO=tf(1,[1 9 27 35])
FTBF=feedback(FTBO,1);
p=rlocus(FTBO,1);
t=[0:0.01:10]';
u=ones(size(t));
figure(1)
rlocus(FTBO)
hold on
plot(real(p), imag(p),'p')
title('lieux des racine')
figure(2)
lsim(FTBF,u,t)
title('reponse à léchelon')
%%
%Q10c
clc
clear
close all
FTBO=tf([1 2],[1 6 12 0])
FTBF=feedback(FTBO,1);
p=rlocus(FTBO,1);
t=[0:0.01:10]';
u=ones(size(t));
figure(1)
rlocus(FTBO)
hold on
plot(real(p), imag(p),'p')
title('lieux des racine')
figure(2)
lsim(FTBF,u,t)
title('reponse à léchelon')
%%
%Q10d
clc
clear
close all
FTBO=tf([1 4 5],[1 9 3 0 0])
FTBF=feedback(FTBO,1);
p=rlocus(FTBO,1);
t=[0:0.01:10]';
u=ones(size(t));
figure(1)
rlocus(FTBO)
hold on
plot(real(p), imag(p),'p')
title('lieux des racine')
figure(2)
lsim(FTBF,u,t)
title('reponse à léchelon')

