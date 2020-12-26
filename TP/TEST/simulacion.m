%% SIMULINK (se ejecuta luego de 'tanques.m') %%
% parametros:
%r   = h3e; % referencia.
dh1 = 1.2;  % perturbacion respecto a equilibrio h1.
dh2 = .9;  % perturbacion respecto a equilibrio h2.
dh3 = .1;  % perturbacion respecto a equilibrio h3.
% run simulink
sys3Tanks_v3 % carga sistema en simulink.