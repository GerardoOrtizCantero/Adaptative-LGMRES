%modificado por jccf agosto 2024
%modificado Gerado Ortiz noviembre 2024

function [x,flag,relressvec,time,ciclos] = ...
Adaptative_LGMRES(A, b, m, k,lL, tol, maxit, xInitial, varargin)

% Adaptative LGMRES algorithm
    %
    %   LGMRES ("Loose GMRES") is a modified implementation of the restarted
    %   Generalized Minimal Residual Error or GMRES(m) [1], performed by
    %   appending 'k' error approximation vectors to the restarting Krylov
    %   subspace, as a way to preserve information from previous
    %   discarted search subspaces from previous iterations of the method.
    %
    %   Augments the standard GMRES approximation space with approximations
    %   to the error from previous restart cycles as in [1].
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relresvec, time] = ...
    %       lgmres(A, b, m, k, tol, maxit, xInitial)
    %
    %
    %   Input Parameters:
    %   -----------------
    %
    %   A:          n-by-n matrix
    %               Left-hand side of the linear system Ax = b.
    %
    %   b:          n-by-1 vector
    %               Right-hand side of the linear system Ax = b.
    %
    %   m:          int
    %               Restart parameter (similar to 'restart' in MATLAB).
    %
    %   k:          int
    %               Number of error approximation vectors to be appended
    %               to the Krylov search subspace. Default is 3, but values
    %               between 1 and 5 are mostly used.
    %   lL:         int
    %               Cantidad de errores de aproximacion utilizadas para 
    %               enriquecer el subespacio de Krylov
    %
    %   tol:        float, optional
    %               Tolerance error threshold for the relative residual norm.
    %               Default is 1e-6.
    %
    %   maxit:      int, optional
    %               Maximum number of outer iterations.
    %
    %   xInitial:   n-by-1 vector, optional
    %               Vector of initial guess. Default is zeros(n, 1).
    %
    %   Output parameters:
    %   ------------------
    %
    %   x:          n-by-1 vector
    %               Approximate solution to the linear system.
    %
    %   flag:       boolean
    %               1 if the algorithm has converged, 0 otherwise.
    %
    %   relressvec: (1 up to maxit)-by-1 vector
    %               Vector of relative residual norms of every outer iteration
    %               (cycles). The last relative residual norm is simply given
    %               by relresvec(end).
    %
    %   mvec:       (1 up to maxit)-by-1 vector
    %               Vector of restart parameter values. In case the
    %               unrestarted algorithm is invoked, mvec = NaN.
    %
    %   time:       scalar
    %               Computational time in seconds.
    %
    %   References:
    %   -----------
    %
    %   [1] Baker, A. H., Jessup, E. R., & Manteuffel, T. (2005). A technique
    %   for accelerating the convergence of restarted GMRES. SIAM Journal on
    %   Matrix Analysis and Applications, 26(4), 962-984.
    %
    %   Copyright:
    %   ----------
    %
    %   This file is part of the KrySBAS MATLAB Toolbox.
    %
    %   Copyright 2023 CC&MA - NIDTec - FP - UNA
    %
    %   KrySBAS is free software: you can redistribute it and/or modify it under
    %   the terms of the GNU General Public License as published by the Free
    %   Software Foundation, either version 3 of the License, or (at your
    %   option) any later version.
    %
    %   KrySBAS is distributed in the hope that it will be useful, but WITHOUT
    %   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    %   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    %   for more details.
    %
    %   You should have received a copy of the GNU General Public License along
    %   with this file.  If not, see <http://www.gnu.org/licenses/>.
    %         
 % ----> Sanity check on the number of input parameters
    if nargin < 2
        error("Too few input parameters. Expected at least A and b.");
    elseif nargin > 8
        error("Too many input parameters.");
    end
 % ----> Sanity checks on matrix A
    % Check whether A is non-empty
    if isempty(A)
        error("Matrix A cannot be empty.");
    end
 % Check whether A is square
    [rowsA, colsA] = size(A);
    if rowsA ~= colsA
        error("Matrix A must be square.");
    end

    n = rowsA;
    clear rowsA colsA;
 % ----> Sanity checks on vector b
    % Check whether b is non-empty
    if isempty(b)
        error("Vector b cannot be empty.");
    end
 % Check whether b is a column vector
    [rowsb, colsb] = size(b);
    if colsb ~= 1
        error("Vector b must be a column vector.");
    end
 % Check whether b has the same number of rows as b
    if rowsb ~= n
        error("Dimension mismatch between matrix A and vector b.");
    end

    clear rowsb colsb;
 % Special sanity checks for LGMRES here

    % ----> Default value and sanityu checks for m
    if (nargin < 3) || isempty(m)
        m = min(n, 10);
    end
 % ----> If m > n, error message is printed
    if m > n
        error("m must satisfy: 1 <= m <= n.");
    end
 % ----> If m == n, built-in unrestarted gmres will be used
    if m == n
        warning("Full GMRES will be used.");
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b);
        time = toc();
        x = gmres_x;
        if gmres_flag == 0
            flag = 1;
        else
            flag = 0;
        end
        relresvec = resvec ./ resvec(1, 1);
        return
    end
 % ----> If m < n AND k == 0, built-in gmres(m) will be used
    if any(m < n) && any(k == 0)
        warning("GMRES(m) will be used.");
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b, m);
        time = toc();
        x = gmres_x;
        if gmres_flag == 0
            flag = 1;
        else
            flag = 0;
        end
        relresvec = resvec ./ resvec(1, 1);
        return
    end
 % ----> Default value and sanity checks for k
    if (nargin < 4) || isempty(k)
        k = 3;
    end
 % Default value and sanity checks for tol
    if (nargin < 5) || isempty(tol)
        tol = 1e-6;
    end

    if tol < eps
        warning("Tolerance is too small and it will be changed to eps.");
        tol = eps;
    elseif tol >= 1
        warning("Tolerance is too large and it will be changed to 1-eps.");
        tol = 1 - eps;
    end
 % ----> Default value for maxit
    if (nargin < 6) || isempty(maxit)
        maxit = min(n, 10);
    end

    % ----> Default value and sanity checks for initial guess xInitial
    if (nargin < 7) || isempty(xInitial)
        xInitial = zeros(n, 1);
    end

    % Check whether xInitial is a column vector
    [rowsxInitial, colsxInitial] = size(xInitial);
    if colsxInitial ~= 1
        error("Initial guess xInitial is not a column vector.");
    end

    % Check whether x0 has the right dimension
    if rowsxInitial ~= n
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        error(msg);
    end
    clear rowsxInitial colsxInitial;

% Algorithm setup
    restart = 1;
    r0 = b - A * xInitial;
    res(1, :) = norm(r0);
    relresvec(1, :) = (norm(r0) / res(1, 1));
    iter(1, :) = restart;
    d=lL;

% Matrix with the history of approximation error vectors////////////////////
    zMat = zeros(n, k);

% while number_of_cycles <=k, we run GMRES(m + k) only
    tic(); % start measuring CPU time

 % Call MATLAB built-in gmres.//////////////////////////////////////////////
    % Ref. [1], pag. 968, recommends GMRES(m + k)
    % if no enough approximation error vectors are stored yet.
    [x, gmres_flag, ~, ~, resvec] = ...
        gmres(A, b, m + k, tol, 1, [], [], xInitial);

% Update residual norm, iterations, and relative residual vector
    res(restart + 1, :) = resvec(end);
    iter(restart + 1, :) = restart + 1;
    relresvec(size(relresvec, 1) + 1, :) = resvec(end) / res(1, 1);

% First approximation error vector
    zMat(:, restart) = x - xInitial;

% gmres uses a flag system. We only care whether the solution has///////////
    % converged or not
    if gmres_flag ~= 0 % if gmres did not converge
        flag = 0;
        xInitial = x;
        restart = restart + 1;
    else
        flag = 1;
        time = toc();
        return
    end

%Preallocating for speed
w=zeros(n,m+d);
z=zeros(n,1);
ij=1; %Para matriz z
minitial=m;
%%%%%%%
logres(1,:)=(norm(r0)/res(1,1));
iter(1,:)=restart;
miteracion(1,1)=minitial;%//////////////////////////////////////////////
mmin=1;
mmax=n-1; %se puede considerar que no tiene cota superior, antes de las 1000 iteraciones  no logra alcanzar mmax con alpha=-3 y delta=5
mstep=1;
alpha0=2;
delta0=0;

%%%%%

while flag==0
  %%%%%%%%%%%%%%
     if iter(size(iter,1),:) ~=1
         
        [miter]=pdrule(m,minitial,mmin,res,iter(size(iter,1),:),mstep, mmax,alpha0, delta0); %cab
        m=miter(1,1);
        minitial=miter(1,2);
    else
        m=minitial;
    end
    miteracion(iter(size(iter,1),:)+1,1)=m;
    beta=norm(r0);
    v(:,1)=r0/beta;
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
        xm=xInitial + Z;
        r0=b-A*xm;
        miteracion(size(miteracion,1)+1,1)=m;
        res(restart+1,:)=norm(r0);
        iter(restart+1,:)=restart+1;
        %logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
        logres(size(logres,1)+1,:)=norm(r0)/res(1,1);

        if logres(size(logres,1)) <tol
            flag=1;
        else
            xInitial=xm;                        %update and restart
            restart=restart+1;
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
                h2=zeros(s+1,s);
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
        xm=xInitial+V*minimizer;
        iter(restart+1,:)=restart+1;
        r0=b-A*xm;
        miteracion(size(miteracion,1)+1,1)=m;
        res(restart+1,:)=norm(r0);
        iter(restart+1,:)=restart+1;
        logres(size(logres,1)+1,:)=norm(r0)/res(1,1);
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
            xInitial=xm;                        %update and restart
            restart=restart+1;
         end


    end %if restarted

end  %while flag

%retorno de variables
time=toc;     %Imprime tiempo de ejecucion
lastcycle=size(logres,1);
%tiempoc= [lastcycle tiempo];
%time= tiempoc;
ciclos= lastcycle;
relressvec=logres;


