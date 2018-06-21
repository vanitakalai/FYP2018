%% Fault tolerant with RPP with binary search theta
%Author: Vanita K

% clearvars ('-except' ,keepVars{:});clc;
clear;clc
sysgen;

%define variables for rpp
r=100;
alpha=0.1;


%% Non fault tolerant control

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
    
    XA=P*A+Y*C;
    
    %half plane
    XA+XA'+ 2*alpha*P <=0;
    %disk
    %
    [-r*P,XA;
        XA', -r*P]<=0;
    %conic sector
    [sin(theta_x)*(XA+XA'),cos(theta_x)*(XA-XA');
        -cos(theta_x)*(XA-XA'),sin(theta_x)*(XA+XA')]<=0;
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

    BRL<=0;
    
    cvx_end
    
    tempNFTC=[tempNFTC;[theta_x,Gamma]];
    
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
        L_NFTC=inv(P)*Y;
        if(theta_U-theta_L<theta_tol)
            theta_final=theta_x;
            break;
        end
        %if tolerance not reached set current theta a upper bound
        theta_U=theta_x;
    end
end

try
theta_min_NFTC=max(abs(abs(angle(eig(A+L_NFTC*C)))-pi));
theta_final_NFTC=theta_final;
catch
    disp('No feasibility found')
    theta_final_NFTC=nan;
end



%% Fault tolerant control

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
    
    for i=1:bny
        
        XA=P*A+Y*Delta(:,:,i)*C;
        
        %half plane
        XA+XA'+ 2*alpha*P <=0;
        %disk
        [-r*P,XA;
            XA', -r*P]<=0;
        %conic sector
        [sin(theta_x)*(XA+XA'),cos(theta_x)*(XA-XA');
            -cos(theta_x)*(XA-XA'),sin(theta_x)*(XA+XA')]<=0
        
    end
    
        
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

    BRL<=0;
    
    
    cvx_end
    
    tempFTC=[tempFTC;[theta_x,Gamma]];
    
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
        L_FTC=inv(P)*Y;%set controller gain for eigenvalues
        if(theta_U-theta_L<theta_tol)
            theta_final=theta_x;
            break;
        end
        %if tolerance not reached set current theta a upper bound
        theta_U=theta_x;
    end
end

try
    theta_min_FTC=max(abs(abs(angle(eig(A+L_FTC*C)))-pi));
    theta_final_FTC=theta_final;
catch
    disp('No feasibility found')
    theta_final_FTC=nan;
end

%% Fast Fault tolerant control

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
    variable S_plane(ny,ny) diagonal;
    variable S_disk(ny,ny) diagonal;
    variable S_conic(2*ny,2*ny) diagonal;
    variable Diag(ny,ny) diagonal;
    G=[zeros(ny), Diag;
        -Diag, zeros(ny)];
    variable Gamma;
    P>=eye(n);
    minimise (Gamma);
         
        %half plane
        [P*A+A'*P+2*alpha*P, Y+C'*S_plane';
            Y'+S_plane*C, -S_plane-S_plane']<=0;
        %disk
        T1_disk=[-r*P, P*A;
            A'*P, -r*P];
        T2_disk=[zeros(n,ny); C'];
        T3_disk=[Y', zeros(ny,n)];
        
        [T1_disk, T3_disk'+T2_disk*S_disk';
            T3_disk+S_disk*T2_disk', -S_disk-S_disk']<=0;
        %conic sector
        T1_conic=[sin(theta_x)*(P*A+A'*P), cos(theta_x)*(P*A-A'*P);
            -cos(theta_x)*(P*A-A'*P), sin(theta_x)*(P*A+A'*P)];
        T2_conic=[sin(theta_x/2)*C', -cos(theta_x/2)*C';
            cos(theta_x/2)*C', sin(theta_x/2)*C'];
        T3_conic=[cos(theta_x/2)*Y', sin(theta_x/2)*Y';
            -sin(theta_x/2)*Y', cos(theta_x/2)*Y'];
        
        [T1_conic, T3_conic'+T2_conic*S_conic'+ T2_conic*G';
            T3_conic+S_conic*T2_conic'+ G*T2_conic', -S_conic-S_conic']<=0;
     
        BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
        (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
        Cz,D,-Gamma*eye(nz)];
    
        BRL<=0;
        
    cvx_end
    
    tempFFTC=[tempFFTC;[theta_x,Gamma]];
    
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
        L_FFTC=inv(P)*Y;%set controller gain for eigenvalues
        if(theta_U-theta_L<theta_tol)
            theta_final=theta_x;
            break;
        end
        %if tolerance not reached set current theta a upper bound
        theta_U=theta_x;
    end
end

try
    theta_min_FFTC=max(abs(abs(angle(eig(A+L_FFTC*C)))-pi));
    theta_final_FFTC=theta_final;
    
catch
    disp('No feasibility found')
    theta_final_FFTC=nan;
end




