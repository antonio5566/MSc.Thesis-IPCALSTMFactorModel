
%% Set estimation parameters and data choice
mkdir Empirical_Result/GB_result
Krange = [2:6];

for K=Krange
    K
    
if K==Krange(1)
    dataname = 'IPCADATA_FNW36_RNKDMN_CON'
    % Load data
    clearvars -except Krange K dataname
    load([ dataname]);
    % als_opt
    als_opt.MaxIterations       = 10000;
    als_opt.Tolerance           = 1e-6;
end

%% Start estimation

disp(['IPCA_empirical_GB starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);

[GammaBeta_initial,s,v]    = svds(X,K);

%% ALS
disp([' ALS started at ' datestr(clock)])

% Numerical ALS procedure: starting from QSVD guess or previous factor guess

GammaBeta      = GammaBeta_initial;
Factor      = s*v';%ones(K,T);%
Lambda     = mean(Factor,2);
%% Fits and R2s

% Fits
XFITS_GB            = nan(L,T);
XFITS_pred_GB       = nan(L,T);
for t=1:T

    XFITS_GB(:,t)     = W(:,:,t)*GammaBeta*Factor(:,t);
    XFITS_pred_GB(:,t)   = W(:,:,t)*GammaBeta*Lambda;
end
%%
Z(:,:,t)*GammaBeta
%%
% R2s
xret(LOC(:)==0)         = nan;
totalsos                = mySOS(xret); 

XR2_total_GB            = 1 - mySOS( X - XFITS_GB  )/mySOS(X);
XR2_pred_GB             = 1 - mySOS( X - XFITS_pred_GB  )/mySOS(X);
   
disp(['  estimation completed ' datestr(clock)])


%% Save results
save(['Empirical_Result/PCA_result/' dataname '_K' num2str(K)] ...
    , 'xret' , 'W' , 'date' , 'LOC' , 'Nts' ...
    , 'Gamma*' , 'Factor*' , 'Lambda*' ...
    , 'X*' ...
     ); 

disp('XR2_total_GB XR2_pred_GB')
disp(num2str([XR2_total_GB XR2_pred_GB]))


end