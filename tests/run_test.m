load('C:\Users\Canlort\Downloads\git\adaptative-lgmres-dev\data\raefsky3.mat');
addpath('C:\Users\Canlort\Downloads\git\adaptative-lgmres-dev\src')
%variables de prueba
A=Problem.A;%A
b=ones(size(A,1),1); %b
m=29;%m
k=5; %k
lL=5; %lL
tol=1e-6; %tol
maxit=150; %maxit
xInitial=zeros(size(A,1),1); %xInitial=x0
varargin=[]; %varargin

%adaptative lgmres
[x1,flag1,relressvec1,time1,cycles] = Adaptative_LGMRES(A, b, m, k, lL, tol, maxit, xInitial);
%pd_gmres
[x2, flag2, relresvec2, mvec, time2] = pd_gmres(A, b, m, [], [], tol, maxit, xInitial,varargin);
%matlab gmres
tic; % Iniciar medición de tiempo
[x3, flag3, relresvecGMRES, iter3, resvec3] = gmres(A, b, [], tol, maxit, [], [], xInitial);
time3 = toc; % Tiempo de ejecución total 
%lgmres
[x4, flag4, relresvec4, time4] = lgmres(A, b, m, k, tol, maxit, xInitial, varargin);
%gmres_e
[x5, flag5, relresvec5, kdvec, time5] = gmres_e(A, b, m, k, tol, maxit, xInitial, varargin);

%grafica
figure;
%Adaptative LGMRES
semilogy(1:length(relressvec1), relressvec1, '-o', 'LineWidth', 2); hold on;
lastNorm1 = relressvec1(end); % Última norma residual
label1 = sprintf('Adaptative LGMRES (Tiempo: %.2f s, Última norma: %.2e)', time1, lastNorm1);

% PD-GMRES
semilogy(1:length(relresvec2), relresvec2, '-s', 'LineWidth', 2); 
lastNorm2 = relresvec2(end); % Última norma residual
label2 = sprintf('PD-GMRES (Tiempo: %.2f s, Última norma: %.2e)', time2, lastNorm2);

% MATLAB GMRES
semilogy(1:length(resvec3), resvec3/norm(b), '-^', 'LineWidth', 2);
lastNorm3 = resvec3(end) / norm(b); % Última norma residual
label3 = sprintf('MATLAB GMRES (Tiempo: %.2f s, Última norma: %.2e)', time3, lastNorm3);

% LGMRES
semilogy(1:length(relresvec4), relresvec4, '-d', 'LineWidth', 2); 
lastNorm4 = relresvec4(end); % Última norma residual
label4 = sprintf('LGMRES (Tiempo: %.2f s, Última norma: %.2e)', time4, lastNorm4);

% GMRES-E
semilogy(1:length(relresvec5), relresvec5, '-x', 'LineWidth', 2); 
lastNorm5 = relresvec5(end); % Última norma residual
label5 = sprintf('GMRES-E (Tiempo: %.2f s, Última norma: %.2e)', time5, lastNorm5);

% Configuración del gráfico
xlabel('Iteraciones','FontSize', 18);
ylabel('Norma Residual (log)','FontSize', 18);
title('Convergencia del Método Adaptative LGMRES y Variantes','FontSize', 18);
grid on;

% Ajustar los límites del eje X para que la gráfica tenga espacio para todas las iteraciones
xlim([1, max([length(relressvec1), length(relresvec2), length(resvec3), length(relresvec4), length(relresvec5)])]);

% Leyenda con tiempos de ejecución y normas directamente
legend({label1, label2, label3, label4, label5}, 'Location', 'northeast', 'FontSize', 20);


