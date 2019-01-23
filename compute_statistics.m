function [re1_s, rate1,var1, FanoFactor, mean_corr]=compute_statistics(s1,N1,Nc, ifCnt, ifFA) 
% s1, s2: spike times of neuron Pop1 & Pop2, respectively, 2xN
% N1, N2: total # of neurons in Pop1 & Pop2, respectively
% Nc: 1x2 # of sampled neurons, e.g. Nc=[500, 500];
% only neuron in the center square [.25, .75]x[.25, .75] are sampled 
% Cbar=[C22, C12, C11], mean correlations 
% C_d: [20x3] correlation for each distance (daxis) 
% COV_d, COVbar: same as C_d, Cbar for covariance 
% rate1, rate2: mean rate (Hz) of the Nc sampled neurons 
% var1, var2: variance of the Nc sampled neurons 

T=floor(min(max(s1(1,:))));

s1=s1(:,s1(1,:)<=T);

N11=round(sqrt(N1));
if size(s1,1)==3
s1(2,:)=(s1(2,:)-1)*N11+s1(3,:); % x=ceil((I)/Nx1), y=mod(I-1,Nx1)+1
end

I1=transpose(unique(s1(2,:)));

I1=I1(I1<=N1);

Ix10=(ceil(I1/N11))/N11;
Iy10=(mod((I1-1),N11)+1)/N11;

I1=I1(Ix10<0.75 & Ix10>0.25 & Iy10<0.75 & Iy10>0.25);

Ic1=randsample(I1,Nc(1));
% compute spike counts using sliding window 
Tw=200; % sliding window size 
Tburn=1000; 
time=0:1:T;

re1=zeros(Nc(1),length(time));

for mm=1:Nc(1)
    re1(mm,:)=hist(s1(1,Ic1(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1(mm)+1/4),time)*1e3;
end

re1_s=imfilter(re1(:,Tburn+1:end),ones(1,Tw)/Tw);re1_s=re1_s(:,Tw/2-1:end-Tw/2);

ind1=mean(re1_s,2)>2;
re1_s=re1_s(ind1,:);

Ix1=(ceil((Ic1(ind1))/N11))/N11;
Iy1=(mod((Ic1(ind1)-1),N11)+1)/N11;
Nc(1)=length(Ix1);

D = pdist2([[Ix1;Ix1],[Iy1;Iy1]],[[Ix1;Ix1],[Iy1;Iy1]],'euclidean');

      

COV=cov([re1_s; re1_s]');
Var=diag(COV);
var1=Var(1:Nc);
rate1=mean(re1_s,2);

%fano factor
FanoFactor = mean(var1./rate1);

R = COV./sqrt(Var*Var'); 

dmax=0.5;
D=D(1:500,1:500);
mean_corr=mean(R(D<=dmax));

end



