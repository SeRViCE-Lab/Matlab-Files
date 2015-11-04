
%%
clear
syms s k1 k2 ka real
A=[-0.275  -.05848; 0.3155  0.05955];
B=[-0.00016;  8.9e-5];
C=[0.742,  2.316];

[Klqr, S, E] = lqr(A,B,diag([1.52,8.11]), 1.86);

Alqr=A-B*Klqr;
eiglqr=eig(Alqr);2

sysc=ss(Alqr,B,C,0);
sysd = c2d(sysc,0.1667/3,'zoh')
eiglqr=eig(sysd.a);


%% design robust system and solve for control gains
Ar = [sysd.a+sysd.b*[k1 k2]    sysd.b*ka;  -sysd.c  0];  
%find d(s) as function of gains
dsk=det(s*eye(3)-Ar);
dsk=collect(dsk);
vpa(dsk,3)

% find desired d(s) 
dsd=(s-eiglqr(1))*(s-eiglqr(1))*(s-.95);
dsd=collect(dsd);
vpa(dsd,3)

%%  solve equations
%solve for k's
Ks=solve( (8.83e-6*k1 - 4.88e-6*k2 - 1.99)==2.93, ...
(4.96e-6*k2 - 8.84e-6*k1 + 4.74e-6*ka + 0.988)==2.86,...
-4.92e-6*ka==-0.929);
Kr= double([Ks.k1, Ks.k2, Ks.ka])