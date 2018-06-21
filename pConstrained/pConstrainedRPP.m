%% p-constrained with Hinfinity method with RPP
%Author: Vanita K
% 
clearvars ('-except' ,keepVars{:}) ;clc;
% 
%generate system param
sysgen;

gamma_RPP=zeros(1,ny+1);
gamma_RPP_relaxed=zeros(1,ny+1);
F_RPP=zeros(ny+1,ny);
F_RPP_relaxed=zeros(ny+1,ny);

%rpp params
r=100;
alpha=1;
theta=(3/4)*pi/2;

%% Fault Tolerant with RPP, p-constrained
%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable Gamma;
    minimize (Gamma);
    P>=eye(n);
    
    %iterate through all fault matrices
    for k=1:bny
        
        % check number of working sensors in each fault matrix>=p
        if int16(trace(Delta(:,:,k)))>=p
            
            AX=P*A+Y*Delta(:,:,k)*C;
            
            %half plane
            AX+AX'+ 2*alpha*P <=0;
            %disk
            [-r*P,AX;
                AX', -r*P]<=0;
            %conic sector
            [sin(theta)*(AX+AX'),cos(theta)*(AX-AX');
                -cos(theta)*(AX-AX'),sin(theta)*(AX+AX')]<=0
        end
    end
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
 
    cvx_end
    
    gamma_RPP(p+1)=Gamma;
    cvx_RPP(p+1)=string(cvx_status);
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    %finds the number of faulty sensors in fault scenarios that lead to instability 
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP,2)
            F_RPP(p+1,i)=sum(x1==i);
        end
    end
    
end


%% Fault Tolerant with relaxed RPP, p-constrained

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable M_disk(ny,ny) symmetric;
    variable M_plane(ny,ny) symmetric;
    variable M_conic(2*ny,2*ny) symmetric;
    variable S_disk(ny,ny) diagonal
    variable S_plane(ny,ny) diagonal
    variable S_conic(2*ny,2*ny) diagonal;
    variable Diag(ny,ny) diagonal;
    variable Gamma;
    minimize (Gamma);
    P>=eye(n);
    M_disk>=0;
    M_plane>=0;
    M_conic>=0;
    G=[zeros(ny), Diag;
    -Diag, zeros(ny)];
    R=(ny-p)*eye(ny)-ones(ny);
    R_conic=2*(ny-p)*eye(2*ny)-ones(2*ny);
    
   %relaxed formulation of iterative method with p
    
   %half plane
   T1_plane=P*A+A'*P+2*alpha*P;
   T2_plane=C';
   T3_plane=Y';
   
   plane_relaxed=[T1_plane+T2_plane*(M_plane.*R)*T2_plane', T3_plane'+T2_plane*(S_plane'-M_plane.*R);
        T3_plane+(S_plane-M_plane.*R)*T2_plane', -S_plane-S_plane'+(M_plane.*R)];
   plane_relaxed<=0;
   %disk
    T1_disk=[-r*P, P*A;
        A'*P, -r*P];
    T2_disk=[zeros(n,ny); C'];
    T3_disk=[Y', zeros(ny,n)];
    
    disk_relaxed=[T1_disk + T2_disk*(M_disk.*R)*T2_disk', T3_disk'+T2_disk*(S_disk'-M_disk.*R); 
        T3_disk+(S_disk-M_disk.*R)*T2_disk', -S_disk-S_disk'+(M_disk.*R)];
    disk_relaxed<=0;
    
   %conic sector
    T1_conic=[sin(theta)*(P*A+A'*P), cos(theta)*(P*A-A'*P);
            -cos(theta)*(P*A-A'*P), sin(theta)*(P*A+A'*P)];
    T2_conic=[sin(theta/2)*C', -cos(theta/2)*C';
            cos(theta/2)*C', sin(theta/2)*C'];
    T3_conic=[cos(theta/2)*Y', sin(theta/2)*Y';
            -sin(theta/2)*Y', cos(theta/2)*Y'];
        
    conic_relaxed=[T1_conic+T2_conic*(M_conic.*R_conic)*T2_conic', T3_conic'+T2_conic*(S_conic'+G'-M_conic.*R_conic);
        T3_conic+(S_conic+G-M_conic.*R_conic)*T2_conic', -S_conic-S_conic'+(M_conic.*R_conic)];
    conic_relaxed<=0;
    
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;

    cvx_end
    
    gamma_RPP_relaxed(p+1)=Gamma;
    cvx_RPP_relaxed(p+1)=string(cvx_status);
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    %finds the number of faulty sensors in fault scenarios that lead to instability 
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_relaxed,2)
            F_RPP_relaxed(p+1,i)=sum(x1==i);
        end
    end
end

