close all;
clear all;
clc;

%Parámetros del sistema:

M = 0.912; %[kg] - Masa del carro
m1 = 2.556e-3; %[kg] - Masa de la barra inferior
m2 = 3.444e-3; %[kg] - Masa de la barra superior
L1 = 0.07378;
L2 = 0.09942;
l1 = L1/2; %[m] 
l2 = L2/2; %[m] 
I1 = 1.159e-6;
I2 = 2.837e-6;
g = 9.81; %m/s2

%algunas constantes que simplifican la notacion a futuro
d1 = m1*l1^2+I1+m2*L1^2+M*L1^2;
d2 = m2*L1*l2+M*L1*L2;
d3 = m2*l2^2+I2+M*L2^2;

f1 = (m1*l1+m2*L1+M*L1)*g;
f2 = (m2*l2+M*L2)*g;

%referencia: como el sistema es autoregulable, la referencia es nula
r=0;

%Ahora planteamos las matrices que componen la ecuación de estados:

%Matrices de la dinámica del sistema: A y B
d =[[d1 d2];[d2 d3]];
delta = det(d);
di = inv(d);
g = [[-f1 0];[0 -f2]];
dg =  -di*g;
h = [1 0]';

A = [[zeros(2,2) eye(2)];[dg zeros(2,2)]];
B = [[zeros(2,1)]; [di*h]];
%Vamos a sensar la posición del carro: consigo controlabilidad

C = [1 0 0 0];
D = 0;

%calculamos las matrices de controlabilidad y observabilidad
aa = poly(A); %coeficientes del pol. característico de A

%Conntrolabilidad:
Wr = [B A*B A^2*B A^3*B]; %no canonica

R_wr = rank(Wr); %es de rango 4 - completo - controlable

Wr2 = [ [1 aa(2) aa(3) aa(4) ];
        [0   1   aa(2) aa(3) ];
        [0   0     1   aa(2) ];
        [0   0     0     1   ]];

 Wr_canonica = inv(Wr2); %matriz canonica de controlabilidad
    
%Observabilidad
Wo = [C; C*A; C*A^2;C*A^3];
R_wo = rank(Wo); %es de rango 4 - completo - observable

Wo2 = [ [  1     0     0     0     ];
        [aa(2)   1     0     0     ];
        [aa(3) aa(2)   1     0     ];
        [aa(4) aa(3) aa(2)   1     ]];
    
Wo_canonica = inv(Wo2); %matriz canonica de observabilidad

% polos de respuesta al escalon
eph = 0.8;  %eph = input('eph: ');
wn = 0.7;    %wn  = input('wn: ');
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
pk = [-sigmad+wd*1i, -sigmad-wd*1i, -sigmad*20, -sigmad*23];
% polos del observador:
pl = -floor(sigmad*100).*[1, 1, 1, 1] + [0, -1, -2, -3];
%
K  = acker(A, B, pk);
Kr = -1/((C/(A-B*K))*B); % ganancia Kr 'feedfoward'.
L  = acker(A', -C', pl)';

%% transferencias planta-controlador
% Transferencia planta:
Ps = tf(ss(A,B,C,D));

% Transferencia controlador:
Ac = A - B*K + L*C;                               
Bc = -L;                                          
Cc = K;
Dc = 0;
Cs = tf(ss(Ac,Bc,Cc,Dc));

% Transferencia lazo cerrado:
k=1;
Ls = k*Cs*Ps;
Ts = minreal(Cs*Ps/(1+Cs*Ps));
