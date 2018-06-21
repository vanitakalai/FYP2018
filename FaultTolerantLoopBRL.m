%% Fault tolerant with Lyapunov Stability
%Author: Vanita K

while(true)
clear; clc;
% clearvars ('-except' ,keepVars{:}) %for loop testing

%generate system param
sysgen;

%% Non fault tolerant control
cvx_begin sdp
cvx_precision high
variable Y(n,ny);
variable P(n,n) symmetric;
variable Gamma;
minimize (Gamma);
P>=eye(n);


BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

BRL<=0;

cvx_end

L_NFTC=inv(P)*Y;
gamma_NFTC=Gamma;
cvx_NFTC=cvx_optval;

%% Fault tolerant control
cvx_begin sdp
cvx_precision high
variable Y(n,ny);
variable P(n,n) symmetric;
variable Gamma;
minimize (Gamma);
P>=eye(n);

for i=1:bny
    %lypanouv stability
    P*A+Y*Delta(:,:,i)*C+(P*A+Y*Delta(:,:,i)*C)'<=0;
end

BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

BRL<=0;

cvx_end

L_FTC=inv(P)*Y; %generate controller gain
gamma_FTC=Gamma;
cvx_FTC=cvx_optval;

%% Fast fault tolerant control 
cvx_begin sdp
cvx_precision high
variable Y(n,ny);
variable P(n,n) symmetric;
variable S(ny,ny) diagonal
variable Gamma;
minimize (Gamma);
P>=eye(n);


%lypanouv stability with relaxation
[P*A+(P*A)', Y+C'*S';
   Y'+S*C, -S-S']<=0


BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

BRL<=0;

cvx_end

L_FFTC=inv(P)*Y; %generate controller gain
gamma_FFTC=Gamma;
cvx_FFTC=cvx_optval;

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
%     if(~isinf(cvx_FTC))
%     EigenValuesFTC(:,i)=eig(A+L_FTC*Delta(:,:,i)*C);
%     end
%     if(~isinf(cvx_NFTC))
%     EigenValuesNFTC(:,i)=eig(A+L_NFTC*Delta(:,:,i)*C);
%     end
%     if(~isinf(cvx_FFTC))
%     EigenValuesFFTC(:,i)=eig(A+L_FFTC*Delta(:,:,i)*C);
%     end
% end
% 
% s_NFTC=nan;
% s_FTC=nan;
% s_FFTC=nan;
% if(~isinf(cvx_NFTC))
% s_NFTC_RPP=nnz(sum(EigenValuesNFTC>=0));
% end
% if(~isinf(cvx_FTC))
% s_FTC_RPP=nnz(sum(EigenValuesFTC>=0));
% end
% if(~isinf(cvx_FFTC))
% s_FFTC_RPP=nnz(sum(EigenValuesFFTC>=0));
% end



