clear all;
close all;
clc

ceros=[];
polos=[1+j,1-j];
k=1;

sys=zpk(ceros,polos,k)
figure
bode(sys)