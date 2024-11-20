
load('C:\Users\Canlort\Downloads\git\krysbas-dev\src\matrices de prueba\memplus.mat');

A=Problem.A;%A
b=ones(size(A,1),1); %b
m=29;%m
k=[]; %k
lL=2; %lL
tol=1e-9; %tol
maxit=500; %maxit
xInitial=zeros(size(A,1),1); %xInitial=x0
varargin=[]; %varargin
%Adaptative Lgmres
[x,flag,relressvec,time,ciclos] = Adaptative_LGMRES(A, b, m, k, lL, tol, maxit, xInitial);

% Graficar la norma residual (NormaResidual) en escala semilogarítmica
figure;
semilogy(1:length(relressvec), relressvec, '-o', 'LineWidth', 2);
xlabel('Iteraciones');
ylabel('Norma Residual (log)');
title('Convergencia del Método LGMRES Adaptativo');
grid on;

str = sprintf('Tiempo de ejecución: %.2f segundos', time);
annotation('textbox', [0.65, 0.8, 0.3, 0.1], 'String', str, ...
           'EdgeColor', 'black', 'BackgroundColor', 'white', ...
           'FontSize', 10, 'HorizontalAlignment', 'center');