%% Set estimation parameters and data choice
clear all
%%
dataname = 'IPCADATA_FNW36_RNKDMN_CON_K5_macrosLASSOFLEX'
K=5
clearvars -except Krange K dataname
load(['Empirical_Result/macrofactor/' dataname]);
%%
Rrange              = [30:36];
collect             = nan(size(GBGD_GB,1),1);
bootsims            = 1000;
[N,L,T] = size(Z);
XFITS_GBGD_null            = nan(size(GammaBeta,1),T);
%% Bootstrap
for R=Rrange
als_opt.MaxIterations       = 10000;
als_opt.Tolerance           = 1e-6;

    nuGBGD_GB = GBGD_GB;
    nuGBGD_GB(R, :) = 0;

    parfor t=1:T
    XFITS_GBGD_null(:,t)        = W(:,:,t)*(nuGBGD_GB*Factor(:,t)+W(:,:,t)*GBGD_GD*mf(:,t)) ;
    end

if bootsims>0
    tic
    
    RESID =  X- XFITS_GBGD;
    boot_GB = nan(1000, size(GBGD_GB,2)+size(GBGD_GD,2));

    rng('shuffle')
    parfor boot=1:bootsims
        % bootstrap indexes
        btix    = [];
        block   = 1;        
        tmp     = 0;
        while tmp<T
            fillin = unidrnd(T) + (0:(block-1));
            fillin = fillin(fillin<=T);
            btix(tmp+(1:length(fillin))) = fillin;
            tmp = tmp+length(fillin);
        end
        btix    = btix(1:T);

        dof     = 5;
        tvar    = dof/(dof-2);

        % Data construction
        % bsxfun calculation among all the factors in matrix 
        X_b     = XFITS_GBGD_null + bsxfun(@times, RESID(:,btix) , (1/sqrt(tvar))*random('t',dof,1,T) );% yes wild
        LOC_1= ~isnan(X_b)  

         
        % Estimation
        GB_Old      = [nuGBGD_GB zeros(L,size(mf,1))];
        F_Old       = Factor;

        tol         = 1;
        iter        = 0;
        while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
            [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X_b,Nts,mf);
            tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
            F_Old   = F_New;
            GB_Old  = GB_New;
            iter    = iter+1;
        end
        
        boot_GB(boot,:) = GB_New(R:R,1:end);
        
    %timing.boot.sims = bootsims;
    %timing.boot.time = toc;
    end

A=[GBGD_GB(R:R,:),GBGD_GD(R:R,:)];
result=mean(sum(boot_GB.^2,2)>mySOS(A))
collect(R,1)=result

end
end

