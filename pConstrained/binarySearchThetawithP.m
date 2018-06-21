%% p-constrained with Binary Search Theta, RPP
%Author: Vanita K

% clear;clc;

clearvars ('-except' ,keepVars{:}) ;clc;
sysgen;

% define variables for rpp
r=100;
alpha=0.1;

theta_final_FTCp=nan(1,ny+1);
theta_final_FFTCp=nan(1,ny+1);
F_FTCp=zeros(ny+1,ny);
F_FFTCp=zeros(ny+1,ny);

%% Fault tolerant control with p


for p=0:ny
    % define theta param
    theta_U=pi;%upper theta bound
    theta_L=0;%lower theta bound
    theta_tol=1e-3;
    
    while(true)
        %set theta at midpoint
        theta_x=(theta_L+theta_U)/2;
        
        cvx_begin sdp
        cvx_precision high
        variable Y(n,ny);
        variable P(n,n) symmetric;
        variable Gamma;
        P>=eye(n);
        minimise (Gamma);
        
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
                [sin(theta_x)*(AX+AX'),cos(theta_x)*(AX-AX');
                    -cos(theta_x)*(AX-AX'),sin(theta_x)*(AX+AX')]<=0
            end
        end
        
        BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
            (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
            Cz,D,-Gamma*eye(nz)];
        
        BRL<=0;
        
        cvx_end
        
        %if no solution found
        if(isinf(cvx_optval))
            %if tolerance reached
            if(theta_U-theta_L<theta_tol/1000)
                %set min theta as upper bound
                theta_final=theta_U;
                %exit loop
                break;
            end
            %if infeasible for pi/2, no feasible solution exists
            if(theta_x>=pi/2)
                break;
            end
            %if tolerance not reached set current theta as lower bound
            theta_L=theta_x;
        else
            L_FTCp=inv(P)*Y;%set controller gain for eigenvalues
            if(theta_U-theta_L<theta_tol)
                theta_final=theta_x;
                break;
            end
            %if tolerance not reached set current theta a upper bound
            theta_U=theta_x;
        end
    end
    try
        theta_final_FTCp(p+1)=theta_final;
        for i=1:bny
            if(~isinf(cvx_optval))
                EigenValues(:,i)=eig(A+L_FTCp*Delta(:,:,i)*C);
            end
        end
        
        if(~isinf(cvx_optval))
            x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
            x1=sum(x'==0);
            for i=1:size(F_FTCp,2)
                F_FTCp(p+1,i)=sum(x1==i);
            end
        end
    catch
        disp('No feasibility found')
    end
    
end

%% Fast Fault tolerant control with p

for p=0:ny
    % define theta param
    theta_U=pi;%upper theta bound
    theta_L=0;%lower theta bound
    theta_tol=1e-3;
    
    while(true)
        %set theta at midpoint
        theta_x=(theta_L+theta_U)/2;
        
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
        P>=eye(n);
        M_disk>=0;
        M_plane>=0;
        M_conic>=0;
        G=[zeros(ny), Diag;
            -Diag, zeros(ny)];
        R=(ny-p)*eye(ny)-ones(ny);
        R_conic=2*(ny-p)*eye(2*ny)-ones(2*ny);
        minimise (Gamma);
        
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
        T1_conic=[sin(theta_x)*(P*A+A'*P), cos(theta_x)*(P*A-A'*P);
            -cos(theta_x)*(P*A-A'*P), sin(theta_x)*(P*A+A'*P)];
        T2_conic=[sin(theta_x/2)*C', -cos(theta_x/2)*C';
            cos(theta_x/2)*C', sin(theta_x/2)*C'];
        T3_conic=[cos(theta_x/2)*Y', sin(theta_x/2)*Y';
            -sin(theta_x/2)*Y', cos(theta_x/2)*Y'];
        
        conic_relaxed=[T1_conic+T2_conic*(M_conic.*R_conic)*T2_conic', T3_conic'+T2_conic*(S_conic'+G'-M_conic.*R_conic);
            T3_conic+(S_conic+G-M_conic.*R_conic)*T2_conic', -S_conic-S_conic'+(M_conic.*R_conic)];
        conic_relaxed<=0;
        
         
        BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
            (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
            Cz,D,-Gamma*eye(nz)];
        
        BRL<=0;
        
        cvx_end
        
        %if no solution found
        if(isinf(cvx_optval))
            %if tolerance reached
            if(theta_U-theta_L<theta_tol/1000)
                %set min theta as upper bound
                theta_final=theta_U;
                %exit loop
                break;
            end
            %if infeasible for pi/2, no feasible solution exists
            if(theta_x>=pi/2)
                break;
            end
            %if tolerance not reached set current theta as lower bound
            theta_L=theta_x;
        else
            L_FFTCp=inv(P)*Y;%set controller gain for eigenvalues
            if(theta_U-theta_L<theta_tol)
                theta_final=theta_x;
                break;
            end
            %if tolerance not reached set current theta a upper bound
            theta_U=theta_x;
        end
    end
    try
        theta_final_FFTCp(p+1)=theta_final;
        for i=1:bny
            if(~isinf(cvx_optval))
                EigenValues(:,i)=eig(A+L_FFTCp*Delta(:,:,i)*C);
            end
        end
        
        if(~isinf(cvx_optval))
            x=dec2bin(find(max(real(EigenValues))>=0)-1)-'0';
            x1=sum(x'==0);
            for i=1:size(F_FFTCp,2)
                F_FFTCp(p+1,i)=sum(x1==i);
            end
        end
    catch
        disp('No feasibility found')
    end
    
end




