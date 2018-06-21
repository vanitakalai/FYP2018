clear; clc;
N=100;%input number of loops
cvx_quiet true %silences CVX output
waitbar(0,'Please wait...');

for iLoop=1:N
    keepVars={'iLoop','N','gamma_str','stability_mat'};
    testCode; %input code to be testes
    gamma_str(iLoop,:)=[gamm_NFTC,gamma_FTC,gamma_FFTC];%variables to be saved
    stability_mat=[stability_mat;[s_NFTC,s_FTC,s_FFTC]];
    waitbar(iLoop / N)
end

 