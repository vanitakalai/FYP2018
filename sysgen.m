%generate system parameters

n=5;
nu=4;
ny=6;
nd=1;
nz=1;
bny=2^ny;


M1=randn(n);
M2=randn(n);
M3=randn(n);
A=M3*(-M1'*M1+(M2-M2'))/M3;
B=randn(n,nu);
C=randn(ny,n);
D=zeros(nz,nd);

Cz=randn(nz,n);
Bd=rand(n,nd);
Dd=rand(ny,nd);

Delta=faultsgen(ny);
EigenValues=zeros(n,bny);
EigenValuesNFTC=zeros(n,bny); %no fault tolerant control
EigenValuesFTC=zeros(n,bny); %with fault tolerant control
EigenValuesFFTC=zeros(n,bny); %with fast fault tolerant control

%generate marginally stable A
M1_marg=randn(n-1); M2_marg=randn(n-1); M3_marg=randn(n-1); 
T=randn(n);
A_marg=M3_marg*(-M1_marg'*M1_marg+(M2_marg-M2_marg'))/M3_marg;
A_marg = [A_marg, zeros(n-1,1); zeros(1, n)];
A_marg = T\A_marg*T;

%generate unstable A
M1_uns=randn(n-1); M2_uns=randn(n-1); M3_uns=randn(n-1); 
T_uns=randn(n);
A_uns=M3_uns*(-M1_uns'*M1_uns+(M2_uns-M2_uns'))/M3_uns;
A_uns = [A_uns, zeros(n-1,1);abs(randn(1,n)) ];
A_uns = T_uns\A_uns*T_uns;

A_mat=[A,A_marg,A_uns];