close all;
clear all;
clc;

%Parámetros del sistema:

M0 = 0.045; %[kg] - Masa del carro
M1 = 0.018; %[kg] - Masa de la barra inferior
M2 = 0.016; %[kg] - Masa de la barra superior
L1 = 0.09; %[m] - Largo barra inferior
L2 = 0.015; %[m] - Largo barra superior
J1 = 0.006; %[kg*m2] - Momento de inercia barra inferior
J2 = 0.006; %[kg*m2] - Momento de inercia barra superior
g = 9.8; %m/s2

%algunas constantes que simplifican la notacion a futuro
n0=M0+M1+M2;
n1=M1*L1+M2*L2;
n2=M1*L1^2+M2*L2^2;
n3=M2*L2;
n4=M2*L1*L2;
n5=M2*L2^2;
n6=M1*L1^2+M2*L1^2;

%referencia: como el sistema es autoregulable, la referencia es nula
r=0;

%Ahora planteamos las matrices que componen la ecuación de estados:

%Matrices de la dinámica del sistema: A y B
  
   A = [[   0      1         0         0         0         0];
        [   0      0      -5.9634      0      -1.3247      0];
        [   0      0         0         1         0         0];
        [   0      0      145.6183     0      -28.4268     0];
        [   0      0         0         0         0         1];
        [   0      0      -47.6149     0       91.2210     0]];
 
B = [0 22.07202 0 -151.371141 0 -56.32414]';

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
coefs = poly([-50 -51 -52 -53 -54 -55]);
aux = coefs-aa;
aux2 = aux(2:7); %tiro los coeficientes a0 y p0
L=inv(Wo)*Wo_canonica*aux2';

