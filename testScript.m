
%% State Feedback Gain
K_NFTC=place(A,B,eig(A+L_NFTC*C)/10);
K_FTC=place(A,B,eig(A+L_FTC*C)/10);
if(exist('L_FFTC','var'))
K_FFTC=place(A,B,eig(A+L_FFTC*C)/10);
end
%% Eigenvalues 
for i=1:bny
    EigenValuesFTC(:,i)=eig(A+L_FTC*Delta(:,:,i)*C);
    EigenValuesNFTC(:,i)=eig(A+L_NFTC*Delta(:,:,i)*C);
    if(exist('L_FFTC','var'))
    EigenValuesFFTC(:,i)=eig(A+L_FFTC*Delta(:,:,i)*C);
    end
end

 
%check if stable
EigenValuesFTC(EigenValuesFTC>=0)
if(exist('L_FFTC','var'))
EigenValuesFFTC(EigenValuesFFTC>=0)
end
EigenValuesNFTC(EigenValuesNFTC>=0)


%plot eigenvalues 
hold on
plot(EigenValuesFTC,'*','Color','b')
plot(EigenValuesNFTC,'o','Color','r')
if(exist('L_FFTC','var'))
plot(EigenValuesFFTC,'+','Color','g')
end
if(exist('r','var'))
plotRPP(theta,r,alpha);
end
hold off

%% Simulation
[~,iNFTC]=max(max(real(EigenValuesNFTC)));
[~,iFTC]=max(max(real(EigenValuesFTC)));
[~,iFFTC]=max(max(real(EigenValuesFFTC)));

try 
%NFTC simulation
L=L_NFTC;
K=K_NFTC;
D=Delta(:,:,iNFTC);
simOut=sim('Model12016a');
figure;
plot(Time,Err)
ylim([-5 5])
title('NFTC time response')
xlabel('Time')
catch 
    disp('Error occured in NFT system due to instability');
end
%FTC simulation

L=L_FTC;
K=K_FTC;
D=Delta(:,:,iFTC);
simOut=sim('Model12016a');
figure;
plot(Time,Err)
title('FTC time response')
xlabel('Time')

if(EigenValuesFFTC(EigenValuesFFTC~=0)~=0)
%FFTC simulation
L=L_FFTC;
K=K_FFTC;
D=Delta(:,:,iFFTC);
simOut=sim('Model12016a');
figure;
plot(Time,Err)
title('FFTC time response')
xlabel('Time')
end

