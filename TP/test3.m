close all;
clear all;
clc;

%Parámetros del sistema:

m0 = 1.5; %[kg] - Masa del carro
m1 = 0.5; %[kg] - Masa de la barra inferior
m2 = 0.75; %[kg] - Masa de la barra superior
l1 = 0.5; %[m] - Largo barra inferior
l2 = 0.75; %[m] - Largo barra superior
g = 9.8; %m/s2

%algunas constantes que simplifican la notacion a futuro
d1 = m0 + m1 + m2;
d2 = (m1/2 + m2) * l1;
d3 = m2 * l2/2;
d4 = (m1/3 + m2) * l1^2;
d5 = m2 * l1 * l2/2;
d6 = m2 * l2^2/3;
f1 = (m1/2 + m2) * l1 * g;
f2 = m2 * l2 *g/2;

%referencia: como el sistema es autoregulable, la referencia es nula
r=0;

%Ahora planteamos las matrices que componen la ecuación de estados:

%Matrices de la dinámica del sistema: A y B
  
   A = [[   0      0         0         1         0         0];
        [   0      0         0         0         1         0];
        [   0      0         0         0         0         1];
        [   0   -7.4920    0.7985      0         0         0];
        [   0   74.9266   -33.7147     0         0         0];
        [   0   -59.9373   52.1208     0         0         0]];
 
B = [0 0 0 -0.6070 1.4984 -0.2839]';

%Vamos a sensar la posición del carro: consigo controlabilidad

C = [1 0 0 0 0 0];
D = 0;

aa = poly(A); %coeficientes del pol. característico de A

%Armamos la matriz de controlabilidad y observabilidad para chequear
Wr = [B A*B A^2*B A^3*B A^4*B A^5*B];
R_wr = rank(Wr); %es de rango 6 - completo - controlable

Wr2 = [ [1 aa(2) aa(3) aa(4) aa(5) aa(6)];
        [0   1   aa(2) aa(3) aa(4) aa(5)];
        [0   0     1   aa(2) aa(3) aa(4)];
        [0   0     0     1   aa(2) aa(3)];
        [0   0     0     0     1   aa(2)];
        [0   0     0     0     0     1  ]];
    
Wr_canonica = inv(Wr2);

Wo = [C; C*A; C*A^2;C*A^3;C*A^4;C*A^5];
R_wo = rank(Wo); %es de rango 6 - completo - observable

%Calculamos K
coefs = poly([-20 -15 -10 -4 -3 -2]);
aux = coefs - aa;
aux2 = aux(2:7);
K = aux2*Wr_canonica*inv(Wr);


Wo2 = [ [  1     0     0     0     0   0 ];
        [aa(2)   1     0     0     0   0 ];
        [aa(3) aa(2)   1     0     0   0 ];
        [aa(4) aa(3) aa(2)   1     0   0 ];
        [aa(5) aa(4) aa(3) aa(2)   1   0 ];
        [aa(6) aa(5) aa(4) aa(3) aa(2) 1]];
    
Wo_canonica = inv(Wo2);

%Calculamos L
coefs = poly([-15 -16 -17 -18 -19 -20]);
aux = coefs-aa;
aux2 = aux(2:7); %tiro los coeficientes a0 y p0
L=inv(Wo)*Wo_canonica*aux2';

%% transferencias planta-controlador
% Transferencia planta:
Ps = tf(ss(A,B,C,D));
% Transferencia controlador:
Ac = A - B*K + L*C;                               % ?? Ao =  A - B*K - L*C;
Bc = -L;                                          % ?? Bo =  L;
Cc = K;                                           % ?? Co = -K;
Dc = 0;
Cs = tf(ss(Ac,Bc,Cc,Dc));
% Transferencia lazo cerrado:
Ls = Cs*Ps;
Ts = minreal(Cs*Ps/(1+Cs*Ps));

