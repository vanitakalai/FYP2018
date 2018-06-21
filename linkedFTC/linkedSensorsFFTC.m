%% Fast Fault tolerant Control for Linked Sensors  (all methods)
%Author: Vanita K

clearvars ('-except' ,keepVars{:}) ;clc;

r=100;
alpha=1;
theta=(3/4)*pi/2;

sysgen;

gamma_BRL_relaxed=zeros(1,ny+1);
gamma_RPP_relaxed=zeros(1,ny+1);
gamma_BRL_linked=zeros(1,ny+1);
gamma_RPP_linked=zeros(1,ny+1);
gamma_BRL_linked_2=zeros(1,ny+1);
gamma_RPP_linked_2=zeros(1,ny+1);
gamma_RPP_BHP=zeros(1,ny+1);

F_BRL_relaxed=zeros(ny+1,ny);
F_RPP_relaxed=zeros(ny+1,ny);
F_BRL_linked=zeros(ny+1,ny);
F_RPP_linked=zeros(ny+1,ny);
F_BRL_linked_2=zeros(ny+1,ny);
F_RPP_linked_2=zeros(ny+1,ny);
F_RPP_BHP=zeros(ny+1,ny);


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
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL_relaxed,2)
            F_BRL_relaxed(p+1,i)=sum(x1==i);
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
    
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_relaxed,2)
            F_RPP_relaxed(p+1,i)=sum(x1==i);
        end
    end
end

%% Fault Tolerant with relaxed H infinity, p-constrained, linked sensors

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable S(ny,ny) diagonal
    variable M11(2,2) symmetric;
    variable M12(2,2) symmetric;
    variable M13(2,2) symmetric;
    variable M22(2,2) symmetric;
    variable M23(2,2) symmetric;
    variable M33(2,2) symmetric;
    variable Gamma;
    minimize (Gamma);
    P>=eye(n);
    M>=0;
    R11=(ny-p)*eye(2)-ones(2);
    R44=-ones(2);
    
    %relaxed formulation of lyapunov stabiltiy for all p
    T_1=P*A+(P*A)';
    T_2=C';
    T_3=Y';
    M_BHP=[M11*R11,M12*R44,M13*R44;
           M12'*R44,M22*R11,M23*R44;
           M13'*R44,M23'*R44,M33*R11];
    
    [T_1 + T_2*(M_BHP)*T_2', T_3' + T_2*(S'-M_BHP);
        T_3 + (S-M_BHP)*T_2', -S-S'+(M_BHP) ]<=0
    
      
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
    
    cvx_end
    
    gamma_BRL_linked(p+1)=Gamma;
    L=inv(P)*Y;
    
    for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL_linked,2)
            F_BRL_linked(p+1,i)=sum(x1==i);
        end
    end
end

%% Fault Tolerant with relaxed RPP, p-constrained, linked sensors

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable M11_plane(2,2) symmetric;
    variable M12_plane(2,2) symmetric;
    variable M13_plane(2,2) symmetric;
    variable M22_plane(2,2) symmetric;
    variable M23_plane(2,2) symmetric;
    variable M33_plane(2,2) symmetric;
    variable M11_disk(2,2) symmetric;
    variable M12_disk(2,2) symmetric;
    variable M13_disk(2,2) symmetric;
    variable M22_disk(2,2) symmetric;
    variable M23_disk(2,2) symmetric;
    variable M33_disk(2,2) symmetric;
    variable M11_conic(4,4) symmetric;
    variable M12_conic(4,4) symmetric;
    variable M13_conic(4,4) symmetric;
    variable M22_conic(4,4) symmetric;
    variable M23_conic(4,4) symmetric;
    variable M33_conic(4,4) symmetric;
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
    R11=(ny-p)*eye(2)-ones(2);
    R44=-ones(2);
    R11_conic=2*(ny-p)*eye(4)-ones(4);
    R44_conic=-ones(4);
    
    M_BHP_plane=[M11_plane*R11,M12_plane*R44,M13_plane*R44;
        M12_plane'*R44,M22_plane*R11,M23_plane*R44;
        M13_plane'*R44,M23_plane'*R44,M33_plane*R11];
    M_BHP_disk=[M11_disk*R11,M12_disk*R44,M13_disk*R44;
        M12_disk'*R44,M22_disk*R11,M23_disk*R44;
        M13_disk'*R44,M23_disk'*R44,M33_disk*R11];
    M_BHP_conic=[M11_conic*R11_conic,M12_conic*R44_conic,M13_conic*R44_conic;
        M12_conic'*R44_conic,M22_conic*R11_conic,M23_conic*R44_conic;
        M13_conic'*R44_conic,M23_conic'*R44_conic,M33_conic*R11_conic];
   %relaxed formulation of iterative method with p
    
   %half plane
   T1_plane=P*A+A'*P+2*alpha*P;
   T2_plane=C';
   T3_plane=Y';
  
   
   plane_relaxed=[T1_plane+T2_plane*(M_BHP_plane)*T2_plane', T3_plane'+T2_plane*(S_plane'-M_BHP_plane);
        T3_plane+(S_plane-M_BHP_plane)*T2_plane', -S_plane-S_plane'+(M_BHP_plane)];
   plane_relaxed<=0;
   %disk
    T1_disk=[-r*P, P*A;
        A'*P, -r*P];
    T2_disk=[zeros(n,ny); C'];
    T3_disk=[Y', zeros(ny,n)];
    
    disk_relaxed=[T1_disk + T2_disk*(M_BHP_disk)*T2_disk', T3_disk'+T2_disk*(S_disk'-M_BHP_disk); 
        T3_disk+(S_disk-M_BHP_disk)*T2_disk', -S_disk-S_disk'+(M_BHP_disk)];
    disk_relaxed<=0;
    
   %conic sector
    T1_conic=[sin(theta)*(P*A+A'*P), cos(theta)*(P*A-A'*P);
            -cos(theta)*(P*A-A'*P), sin(theta)*(P*A+A'*P)];
    T2_conic=[sin(theta/2)*C', -cos(theta/2)*C';
            cos(theta/2)*C', sin(theta/2)*C'];
    T3_conic=[cos(theta/2)*Y', sin(theta/2)*Y';
            -sin(theta/2)*Y', cos(theta/2)*Y'];
        
    conic_relaxed=[T1_conic+T2_conic*(M_BHP_conic)*T2_conic', T3_conic'+T2_conic*(S_conic'+G'-M_BHP_conic);
        T3_conic+(S_conic+G-M_BHP_conic)*T2_conic', -S_conic-S_conic'+(M_BHP_conic)];
    conic_relaxed<=0;
    
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;

    cvx_end
    
    gamma_RPP_linked(p+1)=Gamma;
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_linked,2)
            F_RPP_linked(p+1,i)=sum(x1==i);
        end
    end
end


%% Fault Tolerant with relaxed H infinity, p-constrained, linked sensors

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable S(ny,ny) diagonal
    variable M11(3,3) symmetric;
    variable M12(3,3) symmetric;
    variable M22(3,3) symmetric;
    variable Gamma;
    minimize (Gamma);
    P>=eye(n);
    M>=0;
    R11=(ny-p)*eye(3)-ones(3);
    R44=-ones(3);
    
    %relaxed formulation of lyapunov stabiltiy for all p
    T_1=P*A+(P*A)';
    T_2=C';
    T_3=Y';
    M_BHP=[M11*R11,M12*R44;
           M12'*R44,M22*R11];
    
    [T_1 + T_2*(M_BHP)*T_2', T_3' + T_2*(S'-M_BHP);
        T_3 + (S-M_BHP)*T_2', -S-S'+(M_BHP) ]<=0
    
      
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;
    
    cvx_end
    
    gamma_BRL_linked_2(p+1)=Gamma;
    L=inv(P)*Y;
    
    for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
    end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_BRL_linked_2,2)
            F_BRL_linked_2(p+1,i)=sum(x1==i);
        end
    end
end

%% Fault Tolerant with relaxed RPP, p-constrained, linked sensors

%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable M11_plane(3,3) symmetric;
    variable M12_plane(3,3) symmetric;
    variable M22_plane(3,3) symmetric;
    variable M11_disk(3,3) symmetric;
    variable M12_disk(3,3) symmetric;
    variable M22_disk(3,3) symmetric;
    variable M11_conic(6,6) symmetric;
    variable M12_conic(6,6) symmetric;
    variable M22_conic(6,6) symmetric;
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
    R11=(ny-p)*eye(3)-ones(3);
    R44=-ones(3);
    R11_conic=2*(ny-p)*eye(6)-ones(6);
    R44_conic=-ones(6);
    
    M_BHP_plane=[M11_plane*R11,M12_plane*R44;
        M12_plane'*R44,M22_plane*R11];
    M_BHP_disk=[M11_disk*R11,M12_disk*R44;
        M12_disk'*R44,M22_disk*R11];
    M_BHP_conic=[M11_conic*R11_conic,M12_conic*R44_conic,
        M12_conic'*R44_conic,M22_conic*R11_conic];
   %relaxed formulation of iterative method with p
    
   %half plane
   T1_plane=P*A+A'*P+2*alpha*P;
   T2_plane=C';
   T3_plane=Y';
  
   
   plane_relaxed=[T1_plane+T2_plane*(M_BHP_plane)*T2_plane', T3_plane'+T2_plane*(S_plane'-M_BHP_plane);
        T3_plane+(S_plane-M_BHP_plane)*T2_plane', -S_plane-S_plane'+(M_BHP_plane)];
   plane_relaxed<=0;
   %disk
    T1_disk=[-r*P, P*A;
        A'*P, -r*P];
    T2_disk=[zeros(n,ny); C'];
    T3_disk=[Y', zeros(ny,n)];
    
    disk_relaxed=[T1_disk + T2_disk*(M_BHP_disk)*T2_disk', T3_disk'+T2_disk*(S_disk'-M_BHP_disk); 
        T3_disk+(S_disk-M_BHP_disk)*T2_disk', -S_disk-S_disk'+(M_BHP_disk)];
    disk_relaxed<=0;
    
   %conic sector
    T1_conic=[sin(theta)*(P*A+A'*P), cos(theta)*(P*A-A'*P);
            -cos(theta)*(P*A-A'*P), sin(theta)*(P*A+A'*P)];
    T2_conic=[sin(theta/2)*C', -cos(theta/2)*C';
            cos(theta/2)*C', sin(theta/2)*C'];
    T3_conic=[cos(theta/2)*Y', sin(theta/2)*Y';
            -sin(theta/2)*Y', cos(theta/2)*Y'];
        
    conic_relaxed=[T1_conic+T2_conic*(M_BHP_conic)*T2_conic', T3_conic'+T2_conic*(S_conic'+G'-M_BHP_conic);
        T3_conic+(S_conic+G-M_BHP_conic)*T2_conic', -S_conic-S_conic'+(M_BHP_conic)];
    conic_relaxed<=0;
    
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;

    cvx_end
    
    gamma_RPP_linked_2(p+1)=Gamma;
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_linked_2,2)
            F_RPP_linked_2(p+1,i)=sum(x1==i);
        end
    end
end

%% Fault Tolerant with relaxed RPP, p-constrained, BHP

%generate permutation matrix
idx=[1,7,2,8,3,9,4,10,5,11,6,12];
I=eye(2*ny);
Per=I(idx,:);
%iterate through all p
for p=0:ny
    
    cvx_begin sdp
    cvx_precision high
    variable Y(n,ny);
    variable P(n,n) symmetric;
    variable M_disk(ny,ny) symmetric;
    variable M_plane(ny,ny) symmetric;
    variable M11_conic(4,4) symmetric;
    variable M12_conic(4,4) symmetric;
    variable M13_conic(4,4) symmetric;
    variable M22_conic(4,4) symmetric;
    variable M23_conic(4,4) symmetric;
    variable M33_conic(4,4) symmetric;
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
    R11_conic=2*(ny-p)*eye(4)-ones(4);
    R44_conic=-ones(4);
    M_BHP_conic=[M11_conic*R11_conic,M12_conic*R44_conic,M13_conic*R44_conic;
        M12_conic'*R44_conic,M22_conic*R11_conic,M23_conic*R44_conic;
        M13_conic'*R44_conic,M23_conic'*R44_conic,M33_conic*R11_conic];
    
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
            cos(theta/2)*C', sin(theta/2)*C']*Per';
    T3_conic=Per*[cos(theta/2)*Y', sin(theta/2)*Y';
            -sin(theta/2)*Y', cos(theta/2)*Y'];
                
    conic_relaxed=[T1_conic+T2_conic*(M_conic.*R_conic)*T2_conic', T3_conic'+T2_conic*(S_conic'+G'-M_conic.*R_conic);
        T3_conic+(S_conic+G-M_conic.*R_conic)*T2_conic', -S_conic-S_conic'+(M_conic.*R_conic)];
    conic_relaxed<=0;
    
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
    BRL<=0;

    cvx_end
    
    gamma_RPP_relaxed_BHP(p+1)=Gamma;
    L=inv(P)*Y;
    
       for i=1:bny
        if(~isinf(cvx_optval))
            EigenValues(:,i)=eig(A+L*Delta(:,:,i)*C);
        end
       end
    
    if(~isinf(cvx_optval))
        x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
        x1=sum(x'==0);
        for i=1:size(F_RPP_relaxed_BHP,2)
            F_RPP_relaxed_BHP(p+1,i)=sum(x1==i);
        end
    end
end