function [Lambda]=fa_eigen(spikeCount, M, Nc,Nsample) 

% s1, s2: spike times of neuron Pop1 & Pop2, respectively, 2xN
% N1, N2: total # of neurons in Pop1 & Pop2, respectively
% Nc: 1x2 # of sampled neurons, e.g. Nc=[500, 500];
% only neuron in the center square [.25, .75]x[.25, .75] are sampled 
% Cbar=[C22, C12, C11], mean correlations 
% C_d: [20x3] correlation for each distance (daxis) 
% COV_d, COVbar: same as C_d, Cbar for covariance 
% rate1, rate2: mean rate (Hz) of the Nc sampled neurons 
% var1, var2: variance of the Nc sampled neurons 

addpath('./fa_Yu/');

Lambda=zeros(Nsample,M);  % eigenvalues 
idx=randperm(Nc*Nsample);
for k=1:Nsample
    Y=spikeCount(idx((k-1)*Nc+1:k*Nc),:);
    [estParams, sth] = fastfa(Y, M);
    L=estParams.L;
    LL=L*L';
    [sth,D]=eig(LL);
    la=diag(D);
    %[m,I]=max(la);
    la=sort(la,'descend');
    Lambda(k,:)=la(1:M);
    %eigvector1(:,ss,pid)=V(:,I)*sign(sum(V(:,I)));
    %SIGN(ss,pid)=sign(sum(V(:,I)));
                       
            
end





