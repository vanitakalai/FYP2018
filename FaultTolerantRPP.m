%% Fault tolerant with RPP
%Author: Vanita K

 while(true)
clear;clc;

%testing with loop
% clearvars ('-except' ,keepVars{:})

sysgen;

%% Define variables for rpp
r=100;
alpha=1;
theta=(3/4)*pi/2;

%% Non fault tolerant control
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
   [-r*P,XA;
       XA', -r*P]<=0;
   %conic sector
   [sin(theta)*(XA+XA'),cos(theta)*(XA-XA');
       -cos(theta)*(XA-XA'),sin(theta)*(XA+XA')]<=0
   
   BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

   BRL<=0;

cvx_end

L_NFTC=inv(P)*Y;
Gamma_NFTC=Gamma;
cvx_NFTC_RPP=cvx_optval;

%% Fault tolerant control 
cvx_begin sdp
cvx_precision high
variable Y(n,ny);
variable P(n,n) symmetric;
variable Gamma;
P>=eye(n);
minimise(Gamma);

for i=1:bny
    
   XA=P*A+Y*Delta(:,:,i)*C;
   
   %half plane
   XA+XA'+ 2*alpha*P <=0;
   %disk
   [-r*P,XA;
       XA', -r*P]<=0;
   %conic sector
   [sin(theta)*(XA+XA'),cos(theta)*(XA-XA');
       -cos(theta)*(XA-XA'),sin(theta)*(XA+XA')]<=0
   
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
          (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
              Cz,D,-Gamma*eye(nz)];

    BRL<=0;
    
end
    
   
cvx_end

L_FTC=inv(P)*Y; %generate controller gain
Gamma_FTC=Gamma;
cvx_FTC_RPP=cvx_optval;

%% Fast Fault Tolerant Controller with G
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
minimise(Gamma);

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
    T1_conic=[sin(theta)*(P*A+A'*P), cos(theta)*(P*A-A'*P);
            -cos(theta)*(P*A-A'*P), sin(theta)*(P*A+A'*P)];
    T2_conic=[sin(theta/2)*C', -cos(theta/2)*C';
            cos(theta/2)*C', sin(theta/2)*C'];
    T3_conic=[cos(theta/2)*Y', sin(theta/2)*Y';
            -sin(theta/2)*Y', cos(theta/2)*Y'];
        
    [T1_conic, T3_conic'+T2_conic*S_conic'+ T2_conic*G';
        T3_conic+S_conic*T2_conic'+ G*T2_conic', -S_conic-S_conic']<=0;
    
    BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

    BRL<=0;
   

cvx_end

L_FFTC=inv(P)*Y; %generate controller gain
Gamma_FFTC=Gamma;
cvx_FFTC_RPP=cvx_optval;

% check if all cvx has been solved-only move on until solved
cvx_check=[cvx_NFTC cvx_FTC cvx_FFTC];
if(isempty(cvx_check(isinf(cvx_check)))==1)
    break;
end

end

% testing controllers using simulation model
testScript;

%% Loop Testing Code Extras
% %% Eigenvalues 
% for i=1:bny
%     if(~isinf(cvx_FTC_RPP))
%     EigenValuesFTC(:,i)=eig(A+L_FTC*Delta(:,:,i)*C);
%     end
%     if(~isinf(cvx_NFTC_RPP))
%     EigenValuesNFTC(:,i)=eig(A+L_NFTC*Delta(:,:,i)*C);
%     end
%     if(~isinf(cvx_FFTC_RPP))
%     EigenValuesFFTC(:,i)=eig(A+L_FFTC*Delta(:,:,i)*C);
%     end
% end
% 
% s_NFTC_RPP=nan;
% s_FTC_RPP=nan;
% s_FFTC_RPP=nan;
% if(~isinf(cvx_NFTC_RPP))
% s_NFTC_RPP=nnz(sum(EigenValuesNFTC>=0));
% end
% if(~isinf(cvx_FTC_RPP))
% s_FTC_RPP=nnz(sum(EigenValuesFTC>=0));
% end
% if(~isinf(cvx_FFTC_RPP))
% s_FFTC_RPP=nnz(sum(EigenValuesFFTC>=0));
% end



