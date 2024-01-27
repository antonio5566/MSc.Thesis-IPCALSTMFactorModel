%% Set estimation parameters and data choice
clear all
%mkdir Empirical_Result/macrofactor
%mkdir Empirical_Result/significant
%%
for K=1:6 %change here
    clearvars -except K
macrochoices = {'macrosLASSO'};
for j=1:length(macrochoices)
    macrochoice = macrochoices{j}

if j==1
dataname = 'IPCADATA_FNW36_LASSO_RNKDMN_CON_199001-201405'
% Load data
load([dataname]);
Nts = sum(LOC);
% als_opt
als_opt.MaxIterations       = 5000;
als_opt.Tolerance           = 1e-6;
end

%% macro

csv_in='macro_factor_LASSO.csv'; %2_16_6
macro_factor=importdata(csv_in);
mfactors=macro_factor.data(1:end,1:6)

%%
mfactor_0=mfactors(:,1)
mfactor_1=mfactors(:,2)
mfactor_2=mfactors(:,3)
mfactor_3=mfactors(:,4)
mfactor_4=mfactors(:,5)
mfactor_5=mfactors(:,6)



mdate=macro_factor.textdata(2:end,2);
date_m_2=datetime(mdate,'Format','yyyyMM')
date_m_1=string(date_m_2,'yyyyMM')
date_m =double(date_m_1)
%%

[date,loc1,loc2] = intersect(date_m,date);
%%
switch macrochoice 
    case 'macro_14'
        mf_0 = [ mfactor_1(loc1) mfactor_4(loc1)]';
    case 'macrosLASSO'
        mf = [mfactor_0(loc1) mfactor_1(loc1) mfactor_2(loc1) mfactor_3(loc1) mfactor_4(loc1) mfactor_5(loc1)]';

    case 'macros'
        mf_0 = [mfactor_0(loc1) mfactor_1(loc1) mfactor_2(loc1) mfactor_3(loc1) mfactor_4(loc1) mfactor_5(loc1) mfactor_6(loc1) mfactor_7(loc1) mfactor_8(loc1) mfactor_9(loc1) mfactor_10(loc1) mfactor_11(loc1) mfactor_12(loc1) mfactor_13(loc1) mfactor_14(loc1) mfactor_15(loc1) ]';
        [idx c sumd]=kmeans(mf_0,size(mf_0,1));
        [coeff, score, latent, ~, explained, mu] = pca(mf_0', 'Rows', 'complete');
        cumulative_sum = cumsum(explained);
        threshold = 95;
        loc0 = find(cumulative_sum >= threshold, 1, 'first');
        pca_mf_1 = score(:, 1:loc0);
        mf=pca_mf_1';
end



mf_L = mean(mf,2);
xret = xret(:,loc2);
Q = Q(:,loc2);
X = X(:,loc2);
W = W(:,:,loc2);
LOC = LOC(:,loc2);
Z = Z(:,:,loc2);
Nts = Nts(loc2);
[N,L,T] = size(Z);


%% Start estimation

disp(['IPCA_empirical_GBGDFF starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);

[GammaBeta_QSVD,s,v]    = svds(X,K);

%% ALS for GB-model
disp([' ALS started at ' datestr(clock)])

% Numerical ALS procedure: starting from QSVD guess
tic;
GB_Old      = GammaBeta_QSVD;
F_Old       = s*v';%ones(K,T);%
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
timing.als_xsvd.time = toc;
timing.als_xsvd.iter = iter;
timing.als_xsvd.tols = tols;
disp([' NUM_ALS for GB-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_xsvd.time) ' sec'])

%% ALS for GBGD-model
% Numerical ALS procedure: starting from GB model GammaBeta guess
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
% size(FF,1)=1,3,4,6
GBGD_GD             = GB_New(:,K+1:K+size(mf,1));
GBGD_F              = F_New;
GBGD_L              = mean(F_New,2);
timing.als_gbgd.time = toc;
timing.als_gbgd.iter = iter;
timing.als_gbgd.tols = tols;
disp([' NUM_ALS for GBGD-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_gbgd.time) ' sec'])

%% ALS for GD-model
% Numerical ALS procedure
tic;
GB_Old      = [zeros(L,size(mf,1))];
F_Old       = [];
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
GD_GD             = GB_New;
timing.als_gd.time = toc;
timing.als_gd.iter = iter;
timing.als_gd.tols = tols;
disp([' NUM_ALS for GD-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_gd.time) ' sec'])


%% Fits and R2s

% Fits
RFITS_GB         = nan(N,T);
RFITS_pred_GB    = nan(N,T);
RFITS_GBGD       = nan(N,T);
RFITS_pred_GBGD  = nan(N,T);
RFITS_GD       = nan(N,T);
RFITS_pred_GD  = nan(N,T);
XFITS_GB         = nan(L,T);
XFITS_GD         = nan(L,T);
XFITS_GBGD       = nan(L,T);
XFITS_pred_GB    = nan(L,T);
XFITS_pred_GD    = nan(L,T);
XFITS_pred_GBGD  = nan(L,T);
for t=1:T
    RFITS_GB(:,t)        = Z(:,:,t)*GammaBeta*Factor(:,t) ;
    RFITS_pred_GB(:,t)   = Z(:,:,t)*GammaBeta*Lambda ;
    RFITS_GBGD(:,t)      = Z(:,:,t)*(GBGD_GB*GBGD_F(:,t) + GBGD_GD*mf(:,t)) ;
    RFITS_pred_GBGD(:,t) = Z(:,:,t)*(GBGD_GB*GBGD_L + GBGD_GD*mf_L) ;
    RFITS_GD(:,t)        = Z(:,:,t)*GD_GD*mf(:,t);
    RFITS_pred_GD(:,t)   = Z(:,:,t)*GD_GD*mf_L;
    XFITS_GB(:,t)        = W(:,:,t)*(GammaBeta*Factor(:,t));
    XFITS_pred_GB(:,t)   = W(:,:,t)*(GammaBeta*Lambda);
    XFITS_GD(:,t)        = W(:,:,t)*(GD_GD*mf(:,t));
    XFITS_pred_GD(:,t)   = W(:,:,t)*(GD_GD*mf_L);
    XFITS_GBGD(:,t)      = W(:,:,t)*(GBGD_GD*mf(:,t) + GBGD_GB*GBGD_F(:,t));
    XFITS_pred_GBGD(:,t) = W(:,:,t)*(GBGD_GD*mf_L + GBGD_GB*GBGD_L);
end
QFITS_GB         = GammaBeta*Factor;
QFITS_pred_GB    = repmat(GammaBeta*Lambda,1,T);
QFITS_GBGD       = GBGD_GB*GBGD_F + GBGD_GD*mf;
QFITS_pred_GBGD  = repmat(GBGD_GB*GBGD_L + GBGD_GD*mf_L,1,T);
QFITS_GD         = GD_GD*mf;
QFITS_pred_GD    = repmat(GD_GD*mf_L,1,T);

% R2s
xret(LOC(:)==0)         = nan;
totalsos                = mySOS(xret); 
RR2_total_GB             = 1 - mySOS( xret(LOC(:)) - RFITS_GB(LOC(:))  )/totalsos;
RR2_pred_GB              = 1 - mySOS( xret(LOC(:)) - RFITS_pred_GB(LOC(:))  )/totalsos;
RR2_total_GBGD             = 1 - mySOS( xret(LOC(:)) - RFITS_GBGD(LOC(:))  )/totalsos;
RR2_pred_GBGD              = 1 - mySOS( xret(LOC(:)) - RFITS_pred_GBGD(LOC(:))  )/totalsos;
RR2_total_GD             = 1 - mySOS( xret(LOC(:)) - RFITS_GD(LOC(:))  )/totalsos;
RR2_pred_GD              = 1 - mySOS( xret(LOC(:)) - RFITS_pred_GD(LOC(:))  )/totalsos;

QR2_total_GB             = 1 - mySOS( Q - QFITS_GB  )/mySOS(Q);
QR2_pred_GB              = 1 - mySOS( Q - QFITS_pred_GB  )/mySOS(Q);
QR2_total_GBGD             = 1 - mySOS( Q - QFITS_GBGD  )/mySOS(Q);
QR2_pred_GBGD              = 1 - mySOS( Q - QFITS_pred_GBGD  )/mySOS(Q);
QR2_total_GD             = 1 - mySOS( Q - QFITS_GD  )/mySOS(Q);
QR2_pred_GD              = 1 - mySOS( Q - QFITS_pred_GD  )/mySOS(Q);

XR2_total_GB     = 1 - mySOS( X - XFITS_GB  )/mySOS(X);
XR2_pred_GB      = 1 - mySOS( X - XFITS_pred_GB  )/mySOS(X);
XR2_total_GD     = 1 - mySOS( X - XFITS_GD  )/mySOS(X);
XR2_pred_GD      = 1 - mySOS( X - XFITS_pred_GD  )/mySOS(X);
XR2_total_GBGD   = 1 - mySOS( X - XFITS_GBGD  )/mySOS(X);
XR2_pred_GBGD    = 1 - mySOS( X - XFITS_pred_GBGD  )/mySOS(X);
   
disp(['  estimation completed ' datestr(clock)])


%% Save results
%Empirical_Result/significant
%Empirical_Result/macrofactor
save(['Empirical_Result/macrofactor/' dataname '_K' num2str(K) '_' macrochoice] ...
    , 'xret' , 'W' , 'date' , 'LOC' , 'Nts' ...
    , 'Gamma*' , 'Factor*' , 'Lambda*' ...
    , 'R*' ...
    , 'Q*' ...
    , 'X*' ...
    , 'timing' ...
    , 'mf' ...
    , 'GBGD_*' , 'GD_*','Z','T','mf','-v7.3' ); 

disp('XR2_total_GB XR2_pred_GB XR2_total_GBGD XR2_pred_GBGD')
disp(num2str([XR2_total_GB XR2_pred_GB XR2_total_GBGD XR2_pred_GBGD XR2_total_GD XR2_pred_GD]))
disp('RR2_total_GB RR2_pred_GB RR2_total_GBGD RR2_pred_GBGD')
disp(num2str([RR2_total_GB RR2_pred_GB RR2_total_GBGD RR2_pred_GBGD RR2_total_GD RR2_pred_GD]))

end
end
