%% Fault tolerant Control for Linked Sensors  (all methods)
%Author: Vanita K

clearvars ('-except' ,keepVars{:})

%generate system param
sysgen;
%linked param
l=3;%number of linked blocks
l_array=[2,2,2];% indicates how sensors linked
Delta_L=faultsgenlinked(ny,l,l_array); %linked faults
[~,~,szL]=size(Delta_L);

%linked 2 param
l2=2;%number of linked blocks
l_array2=[3,3];% indicates how sensors linked
Delta_L2=faultsgenlinked(ny,l2,l_array2); %linked faults
[~,~,szL2]=size(Delta_L2);

% generate RPP param
r=100;
alpha=1;
theta=(3/4)*pi/2;

%generate storage vecs
gamma_RPP=zeros(1,ny+1);
gamma_RPP_l=zeros(1,ny+1);
gamma_RPP_l2=zeros(1,ny+1);
gamma_BRL=zeros(1,ny+1);
gamma_BRL_l=zeros(1,ny+1);
gamma_BRL_l2=zeros(1,ny+1);

F_RPP=zeros(ny+1,ny);
F_RPP_l=zeros(ny+1,ny);
F_RPP_l2=zeros(ny+1,ny);
F_BRL=zeros(ny+1,ny);
F_BRL_l=zeros(ny+1,ny);
F_BRL_l2=zeros(ny+1,ny);
%% H infity with Lyapunov Stability, p constrained
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
    L=inv(P)*Y;
    
    for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL,2)
            F_BRL(p+1,i)=sum(x1==i);
        end
    end
end
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
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP,2)
            F_RPP(p+1,i)=sum(x1==i);
        end
    end
    
end
%% H infity with Lyapunov Stability, p constrained, linked sensors(3 blocks)
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
    for k=1:szL
        
        % check number of working sensors in each fault matrix>=p
        if int16(trace(Delta_L(:,:,k)))>=p
            P*A+Y*Delta_L(:,:,k)*C+(P*A+Y*Delta_L(:,:,k)*C)'<=0;
        end
    end
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
    
    cvx_end
    
    gamma_BRL_l(p+1)=Gamma;
    L=inv(P)*Y;
    
    for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL_l,2)
            F_BRL_l(p+1,i)=sum(x1==i);
        end
    end
end
%% Fault Tolerant with RPP, p-constrained, linked sensors(3 blocks)
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
    for k=1:szL
        
        % check number of working sensors in each fault matrix>=p
        if int16(trace(Delta_L(:,:,k)))>=p
            
            AX=P*A+Y*Delta_L(:,:,k)*C;
            
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
    
    gamma_RPP_l(p+1)=Gamma;
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_l,2)
            F_RPP_l(p+1,i)=sum(x1==i);
        end
    end
    
end
%% H infity with Lyapunov Stability, p constrained, linked sensors (2 bocks)
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
    for k=1:szL2
        
        % check number of working sensors in each fault matrix>=p
        if int16(trace(Delta_L2(:,:,k)))>=p
            P*A+Y*Delta_L2(:,:,k)*C+(P*A+Y*Delta_L2(:,:,k)*C)'<=0;
        end
    end
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
    
    cvx_end
    
    gamma_BRL_l2(p+1)=Gamma;
    L=inv(P)*Y;
    
    for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL_l2,2)
            F_BRL_l2(p+1,i)=sum(x1==i);
        end
    end
end
%% Fault Tolerant with RPP, p-constrained, linked sensors (2 blocks)
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
    for k=1:szL2
        
        % check number of working sensors in each fault matrix>=p
        if int16(trace(Delta_L2(:,:,k)))>=p
            
            AX=P*A+Y*Delta_L2(:,:,k)*C;
            
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
    
    gamma_RPP_l2(p+1)=Gamma;
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_l2,2)
            F_RPP_l2(p+1,i)=sum(x1==i);
        end
    end
    
end 