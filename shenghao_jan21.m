%  demo code for two-layer network simulations (related to Fig. 3)
%  calls RF2D3layer.m for network simulation 

clear 

Wtype='broadRec';  taudsyni=8; % (ms)   % Fig. 3Aiv
% Wtype='broadRec';  taudsyni=0.5; % (ms)   % Fig. 3Aii
% Wtype='uniformW'; taudsyni=8; % Fig. 3Aiii
% Wtype='uniformW'; taudsyni=0.5; % Fig. 3Ai
data_folder='data/';   % folder name to store data 
dim='2D'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic


%%%%%%%% set options and parameters to change %%%%%%%%%%%
opt.save=1; % save data 
opt.CompCorr=1; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.Layer1only=1; % 1 for two-layer network, 0 for three-layer network  
opt.loadS1=0;
opt.plotPopR=1; % plot population rate
opt.fixW=0;  % use the same weight matrices for multiple simulations 
%     Wseed1=Wseed1_range(nws); % ransom seed for generating weight matrices
%     Wseed2=Wseed2_range(nws); 

dt=.05; % time step size for integration 
T=2000; % total simulation time (ms) 
filename=strrep(sprintf('%sRF%s2layer_%s_tausyni%.03g_test',...
            data_folder,dim, Wtype, taudsyni),'.','d'),

% parameters to change, default parameters are in RF2D3layer.m
taudsyni13 = {4,6,8,10,12};
taursyni13 = {0.6,0.8,1,1.2,1.4};

for i = 1:5
    for j = 1:5
        taudsyni = taudsyni13{i};
        taursyni = taursyni13{j};
        filename=strrep(sprintf('test_res',...
            data_folder,dim, Wtype, taudsyni,taursyni),'.','d'),
        
        ParamChange={'filename', filename;...
            'dt', dt; 'T', T; 'Nc',Nc;'param(1).taudsyn(3)', taudsyni;'param(1).taursyn(3)', taursyni}; 

        RF2D3layer_shenghao(opt, ParamChange) 
        
    end
end
    


%%%test if 500 if enough to generate statistics%%%
fanos500 = zeros(100,1);
mean_corrs500 = zeros(100,1);

for i = 1:100
    [f, m]=compute_statistics(s1,param(1).Ne,Nc);
    fanos500(i)=f;
    mean_corrs500(i)=m;
end

%%%test if 1500 if enough to generate statistics%%%
fanos1500 = zeros(100,1);
mean_corrs1500 = zeros(100,1);

for i = 1:100
    [f, m]=compute_statistics(s1,param(1).Ne,1000);
    fanos1500(i)=f;
    mean_corrs1500(i)=m;
end

save('fano_corr.mat','fanos500','mean_corrs500','fanos1500','mean_corrs1500')



%%%%%%%%%%%%%%%%%%%%%%%%%jan 22 night
rate5=zeros(100,1);
fanos500 = zeros(100,1);
mean_corrs500 = zeros(100,1);

for i = 1:100
    [re1_s, rate1,var1, FanoFactor, mean_corr] = compute_statistics(s1,N1,Nc);
    rate5(i)=mean(rate1);
    fanos500(i)=FanoFactor;
    mean_corrs500(i)=mean_corr;
end
%%
Lambda_test=fa_eigen(re1_s,5,50,9);
Lambda_percent = Lambda_test./sum(Lambda_test,2);
Lambda_percent_avg = mean(Lambda_percent,1);

Lambda_test_=fa_eigen(re1_s,5,50,9);
Lambda_percent_ = Lambda_test_./sum(Lambda_test_,2);
Lambda_percent_avg_ = mean(Lambda_percent_,1);


%%%test if 3000 if enough to generate statistics%%%
fanos3000 = zeros(100,1);
mean_corrs3000 = zeros(100,1);

for i = 1:100
    [f, m]=compute_statistics(s1,param(1).Ne,3000);
    fanos3000(i)=f;
    mean_corrs3000(i)=m;
end

save('fano_corr.mat','fanos500','mean_corrs500','fanos1000','mean_corrs1000','fanos3000','mean_corrs3000')

%%%%%%%% run simulation %%%%%%%%%%%
RF2D3layer(opt, ParamChange) 

%%%%%% plot correlation vs distance %%%%%%%%%%%%%
load(filename,'C','daxis')
figure
plot(daxis,C,'linewidth',1) 
hold on
colororder=lines(3); 
text(.7, 0.9,'Corr (E2-E2)','unit','n','Horiz','l','color',colororder(1,:))
text(.7, 0.8,'Corr (E2-E1)','unit','n','Horiz','l','color',colororder(2,:))
text(.7, 0.7,'Corr (E1-E1)','unit','n','Horiz','l','color',colororder(3,:))
xlabel('distance (a.u.)')
ylabel('correlation')

%%%%%%% raster movie %%%%%%%%%%%%%%%%
load(filename,'s1')
t1=0; % start time (ms)
t2=1000; % end time (ms)
raster2D_ani(s1,t1,t2,200)

