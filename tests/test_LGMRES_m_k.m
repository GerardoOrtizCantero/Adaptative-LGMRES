%modificado por jccf octubre de 2018

load('C:\Users\Canlort\Downloads\git\krysbas-dev\src\matrices de prueba\memplus.mat');

tic;         %Time Control 
tol=1e-9;
maxit=500;
A=Problem.A;
b=ones(size(A,1),1);
color='m'; Name_Matrix='Problem X';
%%b=ones(size(A,1),1);
x0=zeros(size(A,1),1);
%m=mL;
m=29;
 %d=lL;
lL=2;
d=lL;
%s = m + d;
n=size (A,2);
flag=0;
%restart=1;
r=b-A*x0;
res(1,:)=norm(r);
logres(1,:)= norm(r)/res(1,1);
%iter(1,:)=1;
miteracion(1,1)=m;
%Preallocating for speed
w=zeros(n,m+d);
z=zeros(n,1);
ij=1; %Para matriz z

while flag==0
    beta=norm(r);
    v(:,1)=r/beta; 
    h=zeros(m+1,m);
    if size(logres,1)==1
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
        Z=V*minimizer;
        xm=x0 + Z;
        r=b-A*xm;
        miteracion(size(miteracion,1)+1,1)=m;
        %logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
        logres(size(logres,1)+1,:)=norm(r)/res(1,1);

        if logres(size(logres,1)) <tol 
            flag=1;
        else
            x0=xm;                        %update and restart
        end
        
        %Calculo de z(k)
        z(:,ij)= Z;
        %ij=ij+1;
    else
        if ij<=lL
            d=ij;
            ij=ij+1;
        end
        s=m+d;
        for j=1:s                       %modified gram schmidt--Arnoldi
                if j<=m
                    w(:,j)=A*v(:,j);
                else
                    w(:,j)=A*z(:,d-(j-m-1));
                end
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                s=j;
                h2=zeros(s+1,s);    %VERIFICAR!!!...
                for k=1:s
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else        
            v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
        g=zeros(s+1,1);
        g(1,1)=beta;
        for j=1:s                       %plane rotations (QR decompostion)
            P=eye(s+1);   
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            P(j,j)=cos;
            P(j+1,j+1)=cos;
            P(j,j+1)=sin;
            P(j+1,j)=-sin;
            h=P*h;
            g=P*g;
        end
        R=zeros(s,s);
        G=zeros(s,1);
        V=zeros(n,s);
        for k=1:s
            G(k)=g(k);
            V(:,k)=v(:,k);
            for i=1:s
                R(k,i)=h(k,i);
            end
        end
        for k=m+1:s
            V(:,k)=z(:,d-(k-m-1));
        end
        minimizer=R\G;
        xm=x0+V*minimizer;
        r=b-A*xm;
        miteracion(size(miteracion,1)+1,1)=m;
        logres(size(logres,1)+1,:)=norm(r)/res(1,1);
        %logres(size(logres,1)+1,:)=abs(g(s+1,1)/res(1,1));
        aux=V*minimizer;
        Z=z;
        if size(z,2)<lL
            z(:,size(z,2)+1)=aux;
        else
            for j=2:lL
                z(:,j-1)=Z(:,j);
            end
                z(:,lL)=aux;
        end
        
        

        %if abs(g(s+1,1))/res(1,1) <tol || size(logres,1)==maxit
         if logres(size(logres,1),1) <tol || size(logres,1)==maxit
            flag=1;
         else
            x0=xm;                        %update and restart
         end
 
        
    end %if restarted
   
end  %while flag

tiempo=toc;     %Imprime tiempo de ejecucion
%subplot(1,1,1);
%semilogy(logres,'m')
%hold on
figure;
semilogy(1:length(logres), logres, '-o', 'LineWidth', 2);
xlabel('Iteraciones');
ylabel('Norma Residual (logres)');
title('Convergencia del Método LGMRES ');
grid on;
str = sprintf('Tiempo de ejecución: %.2f segundos', tiempo);
annotation('textbox', [0.65, 0.8, 0.3, 0.1], 'String', str, ...
           'EdgeColor', 'black', 'BackgroundColor', 'white', ...
           'FontSize', 10, 'HorizontalAlignment', 'center');

lastcycle=size(logres,1);
tiempoC= tiempo;
ciclos=lastcycle;
NormaResidual=logres;

