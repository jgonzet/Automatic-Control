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

%% polos de respuesta al escalon
eph = 0.81;  %eph = input('eph: ');
wn = 1.5;    %wn  = input('wn: ');
%
err       = exp(-eph*pi/sqrt(1-eph^2));  % error de sobrepico.
tiempo_s  = 4/(eph*wn);                  % tiempo de establecimiento.
%polo:
sigmad    = eph*wn;            % termino real de polo.
wd        = sqrt(1-eph^2)*wn;  % termino imag de polo.
%
fprintf(['-porcentaje de error de sobrepico (OS): %.2f\n'...
         '-tiempo de establecimiento (Ts): %.2f seg\n'], 100*err, tiempo_s)

%% controlador-observador
% polos del controlador:
pk = [-sigmad+wd*1i, -sigmad-wd*1i, -floor(sigmad*10).*[1, 1, 1, 1] + [0, -1, -2, -3]];
% polos del observador:
pl = -floor(sigmad*100).*[1, 1, 1, 1, 1, 1] + [0, -1, -2, -3, -4, -5];
%
K  = acker(A, B, pk);
Kr = -1/((C/(A-B*K))*B); % ganancia Kr 'feedfoward'.
L  = acker(A', -C', pl)';                           %L(1) = 2.3240793652e6;

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