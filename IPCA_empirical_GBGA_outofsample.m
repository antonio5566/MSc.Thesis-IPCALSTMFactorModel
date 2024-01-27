mkdir Empirical_Result/OOS
%% Set estimation parameters and data choice

% % the following can be commented out if using sbatch
clear all
dataname = 'IPCADATA_FNW36_RNKDMN_CON'
for K=1:6
clearvars -except K dataname


% als_opt
als_opt.MaxIterations       = 5000;
als_opt.Tolerance           = 1e-6;
startindex = 60;


%% Start estimation

disp(['IPCA_empirical_GB starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);

% Load data
load([ dataname]);


bigQ    = Q;
bigX    = X;
bigW    = W;
bigNts  = Nts;

%% Estimation

OOSRFITS_pred_GB    = nan(N,T);
OOSXFITS_pred_GB = nan(L,T);
OOSRealFact         = nan(K,T);
OOSRealTan          = nan(T,1);
OOSRFITS_GB         = nan(N,T);
OOSXFITS_GB      = nan(L,T);
OOSQFITS_GB         = nan(L,T);
OOSQFITS_pred_GB    = nan(L,T);
OOSARBPTF           = nan(1,T);


for t=startindex:T-1
    Q       = bigQ(:,1:t);% this is Q known through t
    X       = bigX(:,1:t);% this is X known through t
    W       = bigW(:,:,1:t);% this is W known through t
    Nts     = bigNts(1:t);% this is Nts known through t

    [GammaBeta_XSVD,s,v]    = svds(X,K);
    % Numerical ALS procedure
    tic;
    if t==startindex
        GB_Old      = [GammaBeta_XSVD zeros(L,1)];
        F_Old       = s*v';%ones(K,T);%
    else
        GB_Old      = [GammaBeta GammaAlpha];
        F_Old       = [Factor s*v(end,:)'];
    end
    tol         = 1;
    iter        = 0;
    while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
        [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts,ones(1,t));
        tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
        F_Old   = F_New;
        GB_Old  = GB_New;
        iter    = iter+1;
    end
    GammaBeta  = GB_New(:,1:end-1);
    GammaAlpha = GB_New(:,end);
    Factor     = F_New;
    Lambda     = mean(F_New,2);
    disp([' NUM_ALS for GB-model for t=' num2str(t) ' done after ' num2str(iter) ' iterations'...
        ' at ' datestr(clock) ' after ' num2str(toc) ' sec'])

    % These are truly out-of-sample -- time t known (because Z is lagged such that Z(:,:,t+1) is actually
    % time t information)
    OOSRFITS_pred_GB(:,t+1)     = Z(:,:,t+1)*GammaBeta*Lambda;
    OOSXFITS_pred_GB(:,t+1)     = bigW(:,:,t+1)*GammaBeta*Lambda;
    OOSQFITS_pred_GB(:,t+1)     = GammaBeta*Lambda;
    
    % These are not truly out-of-sample -- because it uses Q info at t+1
    OOSRealFact(:,t+1)          = ( GammaBeta'*bigW(:,:,t+1)*GammaBeta )\( GammaBeta'*bigX(:,t+1) );
    OOSRealTan(t+1)             = tanptfnext(Factor',OOSRealFact(:,t+1)' );
    
    OOSRFITS_GB(:,t+1)          = Z(:,:,t+1)*GammaBeta*OOSRealFact(:,t+1);
    OOSXFITS_GB(:,t+1)          = bigW(:,:,t+1)*GammaBeta*OOSRealFact(:,t+1);
    OOSQFITS_GB(:,t+1)          = GammaBeta*OOSRealFact(:,t+1);
    
    
    tmp = ( GammaAlpha'/(Z(LOC(:,t),:,t)'*Z(LOC(:,t),:,t))*Z(LOC(:,t),:,t)' );
    tmp2 = xret(LOC(:,t),t+1);
    tmp2(isnan(tmp2))=0;
    OOSARBPTF(t+1)    = tmp*tmp2;
    
end


Q = bigQ;
X = bigX;
W = bigW;
Nts = bigNts;

disp(['  estimation completed ' datestr(clock)])


%% Save results
save(['Empirical_Result/OOS/Results_GBGA_outofsample_' dataname '_K' num2str(K)] ...
    , 'xret' , 'W' , 'Q' , 'X' , 'date' , 'LOC' , 'Nts' ...
    , 'OOS*' ); 


end


































