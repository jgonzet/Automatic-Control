%% plots (se ejecuta luego de 'tanques.m')
fprintf('L(s) = C(s)*P(s):')
Ls
eig(Ls)
figure, nyquist(Ls), grid, title('Nyquist: L(s) = C(s)*P(s)')
figure, bode(Ls),    grid, title('Bode: L(s) = C(s)*P(s)')
