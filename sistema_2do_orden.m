clear all;
close all;
clc

wo=0.1;
zeta=0.2;

sis=tf([wo^2],[1,2*wo*zeta,wo^2]);

step(sis)
figure
pzmap(sis)