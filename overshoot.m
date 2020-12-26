clear all;
close all;
clc;

PO=15;
w=1;

zeta=abs(log(PO/100))/(sqrt(pi^2+log(PO/100)^2));

s=tf('s');
second_order=w^2/(s^2+2*zeta*w*s+w^2)

figure
step(second_order)