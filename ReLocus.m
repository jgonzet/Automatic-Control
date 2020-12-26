clear all, close all, clc

s=tf('s');

%si quiero definirlo con coeficientes del pol. caracteristico:
%num=[1];
%den=[1,7,12];
%t1=tf(num,den);

%si quiero definirlo a partir de ceros y polos:
zeros=[];
poles=[+j,-j];

k=-1;
t2=zpk(zeros,poles,k)


rlocus(t2);