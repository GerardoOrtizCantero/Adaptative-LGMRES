%%modificado por jccf abril 2018
%This is a basic GMRES program with restart

load('C:\Users\Canlort\Downloads\git\krysbas-dev\src\matrices de prueba\memplus.mat');

inicio=cputime;         %Time Control 
tol=1e-9;
A=Problem.A;
b=Problem.b;
%b=[1;2;3];
%b=ones(size(A,1),1);
x0=zeros(size(A,1),1);
%m=m_iter;
m=29;
[s,n]=size (A);
%maxit=itermax;
maxit=500;
flag=0;
%color='r'; Name_Matrix='Problem X';
%  bmax=max(max(A));
%  bmin=min(min(A));
%  x0=bmin*ones(n,1)+(bmax-bmin)*rand(n,1);
if (s~=n)
    error ('Matrix not square');
end

[i,j]=size (b);

if (s~=i)
    error ('Vector does not match size of matrix A');
end

if (j~=1)
    error ('Vector is not a column vector')
end

if (size (b)~=size(x0))
    error('Initial guess not right size');
end

restart=1;
r0=b-A*x0;
res(1,:)=norm(r0);
logres(1,:)=norm(r0)/res(1,1);
%logres(1,:)=norm(r0);
%logres=[];
iter(1,:)=restart;
miteracion(iter(size(iter,1),:),1)=m;
AC=[];
Re=[];
Re(:,size(Re,2)+1)=r0;
while flag==0
    %v=zeros(n,m+1);
    %w=zeros(n,m);
    r=b-A*x0;
    beta=norm(r);
    v(:,1)=r/beta; 
    %h=zeros(m+1,m);

    for j=1:m                       %modified gram schmidt--Arnoldi
        w(:,j)=A*v(:,j);
        for i=1:j
            h(i,j)=w(:,j)'*v(:,i);
            w(:,j)=w(:,j)-h(i,j)*v(:,i);
        end
        h(j+1,j)=norm(w(:,j));
        if h(j+1,j)==0
            m=j;
            h2=zeros(m+1,m);    
            for k=1:m
                h2(:,k)=h(:,k);
            end
            h=h2;
        else        
        v(:,j+1)=w(:,j)/h(j+1,j);
        end
    end
    Hs=h;
    g=zeros(m+1,1);
    g(1,1)=beta;
 
    for j=1:m                       %plane rotations (QR decompostion)
        P=eye(m+1);   
        sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
        cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
        P(j,j)=cos;
        P(j+1,j+1)=cos;
        P(j,j+1)=sin;
        P(j+1,j)=-sin;
        h=P*h;
        g=P*g;
    end
    R=zeros(m,m);
    G=zeros(m,1);
    V=zeros(n,m);
    for k=1:m
        G(k)=g(k);
        V(:,k)=v(:,k);
        for i=1:m
            R(k,i)=h(k,i);
        end
    end
    minimizer=R\G;
    xm=x0+V*minimizer;
    r=b-A*xm;
    Re(:,size(Re,2)+1)=r;
 
    
    %res(restart+1,:)=abs(g(m+1,1));
    iter(restart+1,:)=restart+1;
    miteracion(iter(size(iter,1),:),1)=m;
    %logres(restart+1,:)=abs(g(m+1,1)/res(1,1));
    %logres(restart+1,:)=abs(g(m+1,1));
    logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
    %logres(size(logres,1)+1,:)=abs(g(m+1,1));
    
    %if (abs (g(m+1,1))) <tol  || size(logres,1)==maxit
    if (abs (g(m+1,1)))/res(1,1) <tol  || size(logres,1)==maxit    %empleando última componente de g como residuo
         flag=1;
         residuo= (abs (g(m+1,1)))/res(1,1);
    else
        x0=xm;                        %update and restart
        restart=restart+1;
    end
AC(size(AC,1)+1,:)=acos(logres(size(logres,1),1)/logres(size(logres,1)-1,1))*(180/pi());
end
tiempo=cputime - inicio;     %Imprime tiempo de ejecución
% Graficar la norma residual (NormaResidual) en escala semilogarítmica
figure;
semilogy(1:length(logres), logres, '-b', 'LineWidth', 2);
xlabel('Iteraciones');
ylabel('Norma Residual (log)');
title('Convergencia del Método GMRES');
grid on;

str = sprintf('Tiempo de ejecución: %.2f segundos', tiempo);
annotation('textbox', [0.65, 0.8, 0.3, 0.1], 'String', str, ...
           'EdgeColor', 'black', 'BackgroundColor', 'white', ...
           'FontSize', 10, 'HorizontalAlignment', 'center');

lastcycle=size(logres,1)-1;
tiempoC= tiempo;
ciclos= lastcycle;
NormaResidual=logres;