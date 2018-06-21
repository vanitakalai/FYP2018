%% Fault tolerant with itertive loop
%Author: Vanita K

while(true)

clear;clc;

%generate system param
sysgen;

%non fault tolerant control using pole placement
p=-1:-1:-n;
L_NFTC=-place(A',C',p)';
K_NFTC=place(A,B,p/10);

%fault tolerant control 
cvx_begin sdp
cvx_precision high
variable Y(n,ny);
variable P(n,n) symmetric;

P>=eye(n);
P<=10*eye(n);

for i=1:bny
    %lyapunov stability 
    P*A+Y*Delta(:,:,i)*C+(P*A+Y*Delta(:,:,i)*C)'<=0;
end

cvx_end

if (string(cvx_status)~='Infeasible')
    break;
end

end

L_FTC=inv(P)*Y; %generate controller gain
K_FTC=place(A,B,eig(A+L_FTC*C)/10);


for i=1:bny
    EigenValuesFTC(:,i)=eig(A+L_FTC*Delta(:,:,i)*C);
    EigenValuesNFTC(:,i)=eig(A+L_NFTC*Delta(:,:,i)*C);
end

 
%check if stable
EigenValuesFTC(EigenValuesFTC>=0)
EigenValuesNFTC(EigenValuesNFTC>=0)

%plot eigenvalues 
hold on
plot(EigenValuesFTC,'*','Color','r')
plot(EigenValuesNFTC,'o','Color','b')
hold off

%testing controllers using simulation model
testScript;
