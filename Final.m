clear all;
close all;
clc;

ceros=[];
polos=[];
k=1;

G=zpk(ceros,polos,k)


figure
bode(G)
figure
rlocus(G)
figure
nyqlog(G)