%% Hybrid FTC with Lyapunov Stability (iterative and relaxed)
%Author: Vanita K

clearvars ('-except' ,keepVars{:});clc;
% clear;clc;

%generate system param
sysgen;

gamma_HFTC=zeros(ny-1,1);
gamma_HFFTC=zeros(ny-1,1);
err=zeros(ny-1,2);

for a=1:5 %number of active sensors
%decompose active and passive
any=2^a;
DeltaH=Delta(:,:,end-any+1:end);
Delta_p=eye(ny-a); %passive fault matrix
Delta_a=faultsgen(a); %active fault matrix

C_p=C(1:ny-a,:);
C_a=C(end-a+1:end,:);


%% Hybrid fault tolerant control, Iterative Loop with Lyapunov Stability 
cvx_begin sdp
cvx_precision high
variable Y_a(n,a);
variable Y_p(n,ny-a)
variable P(n,n) symmetric;
variable Gamma;
minimize (Gamma)
P>=eye(n);

Y=[Y_a,Y_p];

for i=1:any
      %lyapunov stability 
     P*A+Y_a*Delta_a(:,:,i)*C_a+Y_p*Delta_p*C_p+(P*A+Y_a*Delta_a(:,:,i)*C_a+Y_p*Delta_p*C_p)'<=0
end

BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

BRL<=0;


cvx_end

 L_p=inv(P)*Y_p;
 L_a=inv(P)*Y_a;
 L_FTC=[L_p,L_a];
 
gamma_HFTC(a)=Gamma;



%% Hybid FFTC
cvx_begin sdp
cvx_precision high
variable Y_a(n,a);
variable Y_p(n,ny-a)
variable P(n,n) symmetric;
variable S(a,a) diagonal
variable Gamma;
minimize (Gamma);
P>=eye(n);

Y=[Y_a,Y_p];

[P*A+(P*A)'+Y_p*C_p+(Y_p*C_p)', Y_a+C_a'*S';
   Y_a'+S*C_a, -S-S']<=0


BRL=[P*A+(P*A)'+Y*C+(Y*C)',P*Bd+Y*Dd,Cz';
    (P*Bd+Y*Dd)',-Gamma*eye(nd),D';
    Cz,D,-Gamma*eye(nz)];

BRL<=0;

cvx_end

gamma_HFFTC(a)=Gamma;

L_p=inv(P)*Y_p;
L_a=inv(P)*Y_a;
L_FFTC=[L_p,L_a];


for i=1:any
    EigenValuesHFFTC(:,i)=eig(A+L_FFTC*DeltaH(:,:,i)*C);
    EigenValuesHFTC(:,i)=eig(A+L_FTC*DeltaH(:,:,i)*C);
end

%% run simulation on model and check stable
err_true=zeros(2^ny,1);
L=L_FTC;
for i=1:2^ny
D1=Delta(:,:,i);
DH1=blkdiag(eye(2),Delta(3:end,3:end,i));
try
simOut=sim('modelext');
if (sum(Err(Err>1E5))>0)
    err_true(i)=1;
end
catch
    err_true(i)=1;
end
end

err_true_1=zeros(2^ny,1);
L=L_FFTC;
for i=1:2^ny
D1=Delta(:,:,i);
DH1=blkdiag(eye(2),Delta(3:end,3:end,i));
try
simOut=sim('modelext');
if (sum(Err(Err>1E5))>0)
    err_true_1(i)=1;
end
catch
    err_true_1(i)=1;
end
end


%A+L(Delta_Hat)C stable
EigenValuesHFTC(EigenValuesHFTC>0)
EigenValuesHFFTC(EigenValuesHFFTC>0)
%check sims stable
err(a,:)=[sum(err_true),sum(err_true_1)];

end
