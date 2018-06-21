%% p-constrained with Hinfinity method with Lyapunov Stability
%Author: Vanita K

clearvars ('-except' ,keepVars{:}) ;clc;

%generate system param
sysgen;

gamma_BRL=zeros(1,ny+1);
gamma_BRL_relaxed=zeros(1,ny+1);
F_BRL=zeros(ny+1,ny);
F_BRL_relaxed=zeros(ny+1,ny);

%% Fault Tolerant with H infinity, p-constrained

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable S(ny,ny) diagonal
    variable Gamma;
    minimize (Gamma);
    P>=eye(n);
    
    %iterate through all fault matrices
    for k=1:bny
        
        % check number of working sensors in each fault matrix>=p
        if int16(trace(Delta(:,:,k)))>=p
            P*A+Y*Delta(:,:,k)*C+(P*A+Y*Delta(:,:,k)*C)'<=0;
        end
    end
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
    
    cvx_end
    
    gamma_BRL(p+1)=Gamma;
    cvx_BRL(p+1)=string(cvx_status);
    L=inv(P)*Y;
    %computes all eigenvalues 
    for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    %finds the number of faulty sensors in fault scenarios that lead to instability 
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL,2)
            F_BRL(p+1,i)=sum(x1==i);
        end
    end
end

%% Fault Tolerant with relaxed H infinity, p-constrained

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable M(ny,ny) symmetric;
    variable S(ny,ny) diagonal
    variable Gamma;
    minimize (Gamma);
    P>=eye(n);
    M>=0;
    R=(ny-p)*eye(ny)-ones(ny);
    
    %relaxed formulation of lyapunov stabiltiy for all p
    T_1=P*A+(P*A)';
    T_2=C';
    T_3=Y';
    
    [T_1 + T_2*(M.*R)*T_2', T_3' + T_2*(S'-M.*R);
        T_3 + (S-M.*R)*T_2', -S-S'+(M.*R) ]<=0
    
      
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
    
    cvx_end
    
    gamma_BRL_relaxed(p+1)=Gamma;
    cvx_BRL_relaxed(p+1)=string(cvx_status);
    L=inv(P)*Y;
       %computes all eigenvalues 
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    %finds the number of faulty sensors in fault scenarios that lead to instability 
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL_relaxed,2)
            F_BRL_relaxed(p+1,i)=sum(x1==i);
        end
    end
end
