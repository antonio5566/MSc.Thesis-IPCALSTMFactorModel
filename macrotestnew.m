clear
%collect             = nan(10,6);
collect             = nan(12,6);
%collect             = nan(19,6);

Krange=[3:6]
for K=Krange
   for R=1:12
dataname = ['IPCADATA_FNW36_RNKDMN_CON']
load (['Empirical_Result/macrofactor/' dataname '_K' num2str(K) '_macrosLASSO_macro'  num2str(R)])

%%
clearvars -except  K R Z X XFITS_GBGD GammaBeta Factor W mf Nts T GBGD_F GBGD_GB GBGD_GD collect

%% Bootstrap
bootsims            = 1000;
[N,L,T]             = size(Z);
XFITS_GBGD_null     = nan(L,T);

als_opt.MaxIterations  = 10000;
als_opt.Tolerance      = 1e-6;
nuGBGD_GD              = zeros(size(GBGD_GD));

%%

%%
boot_GB              = nan(bootsims,L);

if bootsims>0
    tic
    
    RESID =  X- XFITS_GBGD;

    for t=1:T
    XFITS_GBGD_null(:,t)        =  W(:,:,t)*(nuGBGD_GD*mf(:,t) + GBGD_GB*GBGD_F(:,t));
    end
    rng('shuffle')

    for boot=1:bootsims
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
         
        % Estimation
        GB_Old  = [GammaBeta zeros(L,size(mf,1))]
        F_Old   = Factor;

        tol         = 1;
        iter        = 0;
        while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
            [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X_b,Nts,mf);
            tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
            F_Old   = F_New;
            GB_Old  = GB_New;
            iter    = iter+1;
        end

        boot_GB  (boot,:) = GB_New(:,K+1:K+1); 
        %boot_GB_2  (boot,:) = GB_New(:,K+2:K+2); 
        %boot_GB_3  (boot,:) = GB_New(:,K+3:K+3); 
        %boot_GB_4  (boot,:) = GB_New(:,K+4:K+4);
        %boot_GB_5  (boot,:) = GB_New(:,K+5:K+5); 
        %boot_GB_6  (boot,:) = GB_New(:,K+6:K+6); 
        %boot_GB_7  (boot,:) = GB_New(:,K+7:K+7);
        %boot_GB_8  (boot,:) = GB_New(:,K+8:K+8); 
        %boot_GB_9  (boot,:) = GB_New(:,K+9:K+9); 
        %boot_GB_10  (boot,:) = GB_New(:,K+10:K+10);     
        %boot_GB_11  (boot,:) = GB_New(:,K+11:K+11); 
        %boot_GB_12  (boot,:) = GB_New(:,K+12:K+12);
        %boot_GB_13  (boot,:) = GB_New(:,K+13:K+13); 
        %boot_GB_14  (boot,:) = GB_New(:,K+14:K+14); 
        %boot_GB_15  (boot,:) = GB_New(:,K+15:K+15); 
        %boot_GB_16  (boot,:) = GB_New(:,K+16:K+16);
        %boot_GB_17  (boot,:) = GB_New(:,K+17:K+17); 
        %boot_GB_18  (boot,:) = GB_New(:,K+18:K+18); 
        %boot_GB_19  (boot,:) = GB_New(:,K+19:K+19);    


    %timing.boot.sims = bootsims;
    %timing.boot.time = toc;
    end
end


%%
result =mean(sum(boot_GB.^2,2)>mySOS(GBGD_GD(:,1:1)))
%result_2=mean(sum(boot_GB_2.^2,2)>mySOS(GBGD_GD(:,2:2)))
%result_3=mean(sum(boot_GB_3.^2,2)>mySOS(GBGD_GD(:,3:3)))
%result_4=mean(sum(boot_GB_4.^2,2)>mySOS(GBGD_GD(:,4:4)))
%result_5=mean(sum(boot_GB_5.^2,2)>mySOS(GBGD_GD(:,5:5)))
%result_6=mean(sum(boot_GB_6.^2,2)>mySOS(GBGD_GD(:,6:6)))
%result_7=mean(sum(boot_GB_7.^2,2)>mySOS(GBGD_GD(:,7:7)))
%result_8=mean(sum(boot_GB_8.^2,2)>mySOS(GBGD_GD(:,8:8)))
%result_9=mean(sum(boot_GB_9.^2,2)>mySOS(GBGD_GD(:,9:9)))
%result_10=mean(sum(boot_GB_10.^2,2)>mySOS(GBGD_GD(:,10:10)))
%result_11=mean(sum(boot_GB_11.^2,2)>mySOS(GBGD_GD(:,11:11)))
%result_12=mean(sum(boot_GB_12.^2,2)>mySOS(GBGD_GD(:,12:12)))
%result_13=mean(sum(boot_GB_13.^2,2)>mySOS(GBGD_GD(:,13:13)))
%result_14=mean(sum(boot_GB_14.^2,2)>mySOS(GBGD_GD(:,14:14)))
%result_15=mean(sum(boot_GB_15.^2,2)>mySOS(GBGD_GD(:,15:15)))
%result_16=mean(sum(boot_GB_16.^2,2)>mySOS(GBGD_GD(:,16:16)))
%result_17=mean(sum(boot_GB_17.^2,2)>mySOS(GBGD_GD(:,17:17)))
%result_18=mean(sum(boot_GB_18.^2,2)>mySOS(GBGD_GD(:,18:18)))
%result_19=mean(sum(boot_GB_19.^2,2)>mySOS(GBGD_GD(:,19:19)))

collect(R,K:K)=result
%collect(2,K:K)=result_2
%collect(3,K:K)=result_3
%collect(4,K:K)=result_4
%collect(5,K:K)=result_5
%collect(6,K:K)=result_6
%collect(7,K:K)=result_7
%collect(8,K:K)=result_8
%collect(9,K:K)=result_9
%collect(10,K:K)=result_10 
%collect(11,K:K)=result_11
%collect(12,K:K)=result_12 
%collect(13,K:K)=result_13
%collect(14,K:K)=result_14 
%collect(15,K:K)=result_15
%collect(16,K:K)=result_16
%collect(17,K:K)=result_17
%collect(18,K:K)=result_18
%collect(19,K:K)=result_19
    end
end
