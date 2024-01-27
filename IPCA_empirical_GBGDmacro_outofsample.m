% mkdir Empirical_Result/OOS

%% Set estimation parameters and data choice
clear all
%mkdir Empirical_Result/macrofactor
%%
K_range = [2:6];
for K=K_range %change here
    clearvars -except K

macrochoices = 'macrosLASSOFLEX';

dataname = 'IPCADATA_FNW36_RNKDMN_CON_KEEPSMA22LL'
% Load data
load([dataname]);
Nts = sum(LOC);

% als_opt
als_opt.MaxIterations       = 5000;
als_opt.Tolerance           = 1e-6;
startindex = 60;


%% macro
csv_in='macro_factor_LASSO_3.csv'; %2_16_6
macro_factor=importdata(csv_in);
mfactors=macro_factor.data(1:end,1:4)
%%
csv_in_1='macro_factor_LASSO_6.csv'
macro_factor_1=importdata(csv_in_1);
mfactors_1=macro_factor_1.data(1:end,1:4)
%%
csv_in_2='macro_factor_LASSO_12.csv'
macro_factor_2=importdata(csv_in_2);
mfactors_2=macro_factor_2.data(1:end,1:4)
%%
%%
mfactor_0=mfactors(:,1)
mfactor_1=mfactors(:,2)
mfactor_2=mfactors(:,3)
mfactor_3=mfactors(:,4)
mfactor_4=mfactors_1(:,1)
mfactor_5=mfactors_1(:,2)
mfactor_6=mfactors_1(:,3)
mfactor_7=mfactors_1(:,4)
mfactor_8=mfactors_2(:,1)
mfactor_9=mfactors_2(:,2)
mfactor_10=mfactors_2(:,3)
mfactor_11=mfactors_2(:,4)

mdate=macro_factor_2.textdata(2:end,2);
date_m_2=datetime(mdate,'Format','yyyyMM')
date_m_1=string(date_m_2,'yyyyMM')
date_m =double(date_m_1)
%%
[date,loc1,loc2] = intersect(date_m,date);

        if K==1
        bigmf= [ ]'   
        elseif K==2
        bigmf= [mfactor_4(loc1)]';
        elseif K==3
        bigmf = [mfactor_11(loc1)]';
        elseif K==4
        bigmf = [mfactor_11(loc1)]';
        elseif K==5
        bigmf = [mfactor_4(loc1)]';
        elseif K==6
        bigmf = [mfactor_7(loc1)]';
        end

xret = xret(:,loc2);
Q = Q(:,loc2);
X = X(:,loc2);
W = W(:,:,loc2);
LOC = LOC(:,loc2);
Z = Z(:,:,loc2);
Nts = Nts(loc2);
[N,L,T] = size(Z);


%% als_opt


%% Start estimation

disp(['IPCA_empirical_GBGDFF starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);
bigX    = X;
bigW    = W;
bigNts  = Nts;


OOSRFITS_pred_GBGD    = nan(N,T);
OOSXFITS_pred_GBGD = nan(L,T);
OOSRealFact         = nan(K,T);
OOSRealFact_1       = nan(K+size(bigmf,1),T);
OOSRealTan          = nan(T,1);
OOSRealTan_1        = nan(T,1);
OOSRFITS_GBGD         = nan(N,T);
OOSXFITS_GBGD      = nan(L,T);

for t=startindex:T-1
    
    X       = bigX(:,1:t);% this is X known through t
    W       = bigW(:,:,1:t);% this is W known through t
    Nts     = bigNts(1:t);%
    mf      = bigmf(:,1:t);
    mf_L     = mean(mf,2);

[GammaBeta_XSVD,s,v]    = svds(X,K);

%% ALS for GB-model
disp([' ALS started at ' datestr(clock)])

% Numerical ALS procedure: starting from QSVD guess
tic;
if t==startindex
        GB_Old      = GammaBeta_XSVD;
        F_Old       = s*v';%ones(K,T);%

else
        GB_Old      = GammaBeta;
        F_Old       = [Factor s*v(end,:)'];
end

tol         = 1;
iter        = 0;
tols        = nan(500,1);
while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
    [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts);
    tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
    F_Old   = F_New;
    GB_Old  = GB_New;
    iter    = iter+1;
    tols    = [tols(2:end);tol];
end
GammaBeta  = GB_New;
Factor     = F_New;
Lambda     = mean(F_New,2);

%timing.als_xsvd.time = toc;
%timing.als_xsvd.iter = iter;
%timing.als_xsvd.tols = tols;
%disp([' NUM_ALS for GB-model done after ' num2str(iter) ' iterations'...
    %' at ' datestr(clock) ' after ' num2str(timing.als_xsvd.time) ' sec'])

%% ALS for GBGD-model
% Numerical ALS procedure: starting from GB model GammaBeta guess
tic;
GB_Old      = [GammaBeta zeros(L,size(mf,1))];
F_Old       = Factor;
tol         = 1;
iter        = 0;
tols        = nan(500,1);


while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
    [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts,mf);
    tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
    F_Old   = F_New;
    GB_Old  = GB_New;
    iter    = iter+1;
    tols    = [tols(2:end);tol];
end
GBGD_GB             = GB_New(:,1:K);
GBGD_GD             = GB_New(:,K+1:K+size(mf,1));
GBGD_F              = F_New;
GBGD_L              = mean(F_New,2);
F= [GBGD_F',mf']  

timing.als_gbgd.time = toc;
timing.als_gbgd.iter = iter;
timing.als_gbgd.tols = tols;
disp([' NUM_ALS for GBGD-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_gbgd.time) ' sec'])

    
    OOSRFITS_pred_GBGD(:,t+1)     = Z(:,:,t+1)*(GBGD_GB*GBGD_L + GBGD_GD*mf_L) 
    OOSXFITS_pred_GBGD(:,t+1)     = bigW(:,:,t+1)*(GBGD_GD*mf_L + GBGD_GB*GBGD_L);
    
    % These are not truly out-of-sample -- because it uses Q info at t+1
    
    OOSRealFact(:,t+1)          = ( GBGD_GB'*bigW(:,:,t+1)*GBGD_GB)\( GBGD_GB'*bigX(:,t+1) )
    OOSRealTan(t+1)             = tanptfnext(F,[OOSRealFact(:,t+1);bigmf(:,t+1)]' )

    OOSRFITS_GBGD(:,t+1)          = Z(:,:,t+1)*(GBGD_GB*OOSRealFact(:,t+1)+ GBGD_GD*bigmf(:,t+1));
    OOSXFITS_GBGD(:,t+1)          = bigW(:,:,t+1)*(GBGD_GB* OOSRealFact(:,t+1) + GBGD_GD* bigmf(:,t+1));

end
       
disp(['  estimation completed ' datestr(clock)])

X = bigX;
W = bigW;
Nts = bigNts;


%% Save results
macrochoices = 'macro_sig'
save(['Empirical_Result/OOS/Results_GBGDmacro_outofsample_' dataname '_' macrochoices '_K' num2str(K) ] ...
    , 'xret' , 'W' , 'date' , 'LOC' , 'Nts' ...
    , 'Q*' ...
    , 'X*' ...
    , 'OOS*' ); 


end

%%
disp(size(GBGD_F));
disp(size(mf));
%%
