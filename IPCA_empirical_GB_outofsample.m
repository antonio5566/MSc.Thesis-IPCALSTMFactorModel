
%% Set estimation parameters and data choice

% the following can be commented out if using sbatch
clear
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
load([dataname]);


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


for t=startindex:T-1
    X       = bigX(:,1:t);% this is X known through t
    W       = bigW(:,:,1:t);% this is W known through t
    Nts     = bigNts(1:t);% this is Nts known through t

    [GammaBeta_XSVD,s,v]    = svds(X,K);
    % Numerical ALS procedure
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
    while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
        [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts);
        tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
        F_Old   = F_New;
        GB_Old  = GB_New;
        iter    = iter+1;
    end
    GammaBeta  = GB_New;
    Factor     = F_New;
    Lambda     = mean(F_New,2);
    disp([' NUM_ALS for GB-model for t=' num2str(t) ' done after ' num2str(iter) ' iterations'...
        ' at ' datestr(clock) ' after ' num2str(toc) ' sec'])

    % These are truly out-of-sample -- time t known (because Z is lagged such that Z(:,:,t+1) is actually
    % time t information)
    OOSRFITS_pred_GB(:,t+1)     = Z(:,:,t+1)*GammaBeta*Lambda;
    OOSXFITS_pred_GB(:,t+1)     = bigW(:,:,t+1)*GammaBeta*Lambda;
    
    % These are not truly out-of-sample -- because it uses Q info at t+1
    OOSRealFact(:,t+1)          = ( GammaBeta'*bigW(:,:,t+1)*GammaBeta )\( GammaBeta'*bigX(:,t+1) );
    OOSRealTan(t+1)             = tanptfnext(Factor',OOSRealFact(:,t+1)' );
    
    OOSRFITS_GB(:,t+1)          = Z(:,:,t+1)*GammaBeta*OOSRealFact(:,t+1);
    OOSXFITS_GB(:,t+1)          = bigW(:,:,t+1)*GammaBeta*OOSRealFact(:,t+1);
end


X = bigX;
W = bigW;
Nts = bigNts;

disp(['  estimation completed ' datestr(clock)])


%% Save results
save(['Empirical_Result/OOS/' 'Results_GB_outofsample_' dataname '_K' num2str(K)] ...
    , 'xret' , 'W' , 'X' ,'Z','Q', 'date' , 'LOC' , 'Nts','GammaBeta','Factor' ...
    , 'OOS*' ,'-v7.3'); 


end

%%
xret(LOC(:)==0)         = nan;
totalsos                = mySOS(xret); 
OOSRFITS_GB_1= 1 - mySOS(xret(LOC(:)) - OOSRFITS_GB(LOC(:)))/totalsos;
OOSRFITS_pred_GB_1= 1 - mySOS(xret(LOC(:)) - OOSRFITS_pred_GB(LOC(:)))/totalsos;
OOSXRFITS_GB_1= 1 - mySOS(X- OOSXFITS_GB)/mySOS(X);
OOSXRFITS_pred_GB_1= 1 - mySOS(X - OOSXFITS_pred_GB)/mySOS(X);




%%
























