% Script complementario del ej. 2.40 NISE 2015/ 2.37 NISE 2010

clear all;
clc;

syms s complex 
% OJO, NO ES LA "s" del toolbox de control s=tf('s')

% Esta es la parte A), resolución en transferencia:

syms JF1Eq JF2Eq K Deq real

Psi = [JF1Eq*s^2+K , -K;...
        -K, JF2Eq*s^2+Deq*s+K];
 
TF_matricial=inv(Psi)

% Transferencia final Theta3/Teq(s)
Theta3SobreTeq=TF_matricial(2,1) % Es el elemento 2,1.


% ESPACIO DE ESTADOS---------------------

Kmonio=K/JF1Eq;
Ksombrero=K/JF2Eq;
Dsombrero=Deq/JF2Eq;

A=[0 0 1 0;...
    0 0 0 1;...
    -Kmonio Kmonio 0 0;...
    Ksombrero -Ksombrero 0 -Dsombrero];
B=[0;0;1/JF1Eq;0];
C=[0 1 0 0];
D=0;

% Transferencia final Theta3/Teq(s), versión 2.
% Notar que esta se está sacando en base a las ecuaciones en espacio de
% estados.

Theta3SobreTeqV2=C*inv(s*eye(4)-A)*B




