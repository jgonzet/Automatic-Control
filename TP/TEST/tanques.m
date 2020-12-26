%% Script - Sistema de 3er orden de TANQUES ACOPLADOS %%
clear all
close all

%% parametros iniciales
S   = .0171; % (m^2)        seccion de area de tanque.
S1  = 5e-5;  % (m^2)        seccion de area de tubo 1.
S2  = 5e-5;  % (m^2)        seccion de area de tubo 2.
S3  = 5e-5;  % (m^2)        seccion de area de tubo 3.
az1 = .511;  %              coeficiente flujo de salida del tubo 1.
az2 = .528;  %              coeficiente flujo de salida del tubo 2.
az3 = .731;  %              coeficiente flujo de salida del tubo 3.
g   = 9.81;  % (m/s^2)      constante de aceleracion gravitacional.
Q2e = .674e-4; % (m^3/s)    constante de caudal de entrada 2.

%% punto de equilibrio
h1e = .6; % (m)             altura 1 en equilibrio.
h2e = .5; % (m)             altura 2 en equilibrio.
h3e = .406; % (m)           altura 3 en equilibrio.
Q1e = .358e-4; % (m^3/s)    caudal de entrada 1 en equilibrio.

%% espacio de estados
a = S1*az1*g/(S*sqrt(2*g*(h1e-h2e)));
b = S2*az2*g/(S*sqrt(2*g*(h2e-h3e)));
c = S3*az3*g/(S*sqrt(2*g*h3e));
% Matrices de Espacio de Estados:
A = [-a, a, 0; a, -a-b, b; 0, b, -b-c];
B = [1/S; 0; 0];
C = [0 0 1];
D = 0;

%% polos de respuesta al escalon
eph = 0.81; % eph = input('eph: ');
wn = 0.2;   % wn  = input('wn: ');
%
err = exp(-eph*pi/sqrt(1-eph^2));  % error de sobrepico.
ts  = 4/(eph*wn);                  % tiempo de establecimiento.
%polo:
sigmad    = eph*wn;            % termino real de polo.
wd        = sqrt(1-eph^2)*wn;  % termino imag de polo.
%
fprintf(['-porcentaje de error de sobrepico (OS): %.2f %c\n'...
         '-tiempo de establecimiento (Ts): %.2f seg\n'], 100*err, '%', ts)

%% controlador-observador
% polos del controlador:
pk = [-sigmad+wd*1i, -sigmad-wd*1i, -sigmad*10];
% polos del observador:
pl = -floor(sigmad*100).*[1, 1, 1] + [0, -1, -2];
%
K  = acker(A, B, pk);
Kr = -1/((C/(A-B*K))*B); % ganancia Kr 'feedfoward'.
L  = acker(A', -C', pl)';                          %L(1) = 2.3240793652e6;

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
