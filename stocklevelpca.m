clc;
clear all

mkdir GB_result/PCA_result

load("IPCADATA_FNW36_RNKDMN_CON.mat");
dataname = 'Results_PCA_RNKDMN_CON'

%% Setting 


%% Start Estimation

disp(['PCA_Regressions starting at ' datestr(clock) ': data=' dataname]);

%% Declare the switch

[date,loc1,loc2] = intersect(date,date);

      RET     =xret(:,loc2);
      Rlogic  = true;   


%% PCA
disp([' PCA started at ' datestr(clock)])

X_norm=zscore(X)
%%
% managed portfolio PCA

[coeff_1, score_1, latent_1,explained1] = pca(X','algorithm','als','NumComponents',1);
[coeff_3, score_3, latent_3,explained3] = pca(X','algorithm','als','NumComponents',3);
[coeff_4, score_4, latent_4,explained4] = pca(X','algorithm','als','NumComponents',4);
[coeff_5, score_5, latent_5,explained5] = pca(X','algorithm','als','NumComponents',5);
[coeff_6, score_6, latent_6,explained6] = pca(X','algorithm','als','NumComponents',6);




%%
n=1:size(RET,1)
%%

FITS_PCA1 = nan(size(RET));
FITS_PCA3 = nan(size(RET));
FITS_PCA4 = nan(size(RET));
FITS_PCA5 = nan(size(RET));
FITS_PCA6 = nan(size(RET));

FITS_cond_PCA1 = nan(size(RET));
FITS_cond_PCA3 = nan(size(RET));
FITS_cond_PCA4 = nan(size(RET));
FITS_cond_PCA5 = nan(size(RET));
FITS_cond_PCA6 = nan(size(RET));

%%
X_reduced_1_mean = mean(X_reduced_1,1);
X_reduced_3_mean  = repmat(mean(X_reduced_3),size(RET,2),1);
X_reduced_4_mean  = repmat(mean(X_reduced_4),size(RET,2),1);
X_reduced_5_mean = repmat(mean(X_reduced_5),size(RET,2),1);
X_reduced_6_mean  = repmat(mean(X_reduced_6),size(RET,2),1);
%%

for n=1:size(RET,1)
               
beta1 = regress(RET(n,:)',X_reduced_1);
beta3 = regress(RET(n,:)',X_reduced_3);
beta4 = regress(RET(n,:)',X_reduced_4);
beta5 = regress(RET(n,:)',X_reduced_5);
beta6 = regress(RET(n,:)',X_reduced_6);

%%

FITS_PCA1(n,:) = X_reduced_1*beta1;
FITS_PCA3(n,:) = X_reduced_3*beta3;
FITS_PCA4(n,:) = X_reduced_4*beta4;
FITS_PCA5(n,:) = X_reduced_5*beta5;
FITS_PCA6(n,:) = X_reduced_6*beta6;
%%

FF= latent_1*beta1;
FITS_cond_PCA3(n,:) = X_reduced_3_mean*beta3;
FITS_cond_PCA4(n,:) = X_reduced_4_mean*beta4;
FITS_cond_PCA5(n,:) = X_reduced_5_mean*beta5;
FITS_cond_PCA6(n,:) = X_reduced_6_mean*beta6;


end
 
%%   
XFITS_GB            = nan(L,T);
for t=1:T

    XFITS_GB(:,t)     = beta1*latent_1;

end
%%
%%
TR2_K1=1-(nansum([RET-FITS_PCA1].^2,"All")/nansum(RET.^2,"All"))
TR2_K3=1-(nansum([RET-FITS_PCA3].^2,"All")/nansum(RET.^2,"All"))
TR2_K4=1-(nansum([RET-FITS_PCA4].^2,"All")/nansum(RET.^2,"All"))
TR2_K5=1-(nansum([RET-FITS_PCA5].^2,"All")/nansum(RET.^2,"All"))
TR2_K6=1-(nansum([RET-FITS_PCA6].^2,"All")/nansum(RET.^2,"All"))
PR2_K1=1-(nansum([RET-FF].^2,"All")/nansum(RET.^2,"All"))
PR2_K3=1-(nansum([RET-FITS_cond_PCA3].^2,"All")/nansum(RET.^2,"All"))
PR2_K4=1-(nansum([RET-FITS_cond_PCA4].^2,"All")/nansum(RET.^2,"All"))
PR2_K5=1-(nansum([RET-FITS_cond_PCA5].^2,"All")/nansum(RET.^2,"All"))
PR2_K6=1-(nansum([RET-FITS_cond_PCA6].^2,"All")/nansum(RET.^2,"All"))

%%
% save(['GB_result/PCA_result/Results_PCA_' XorR '_' dataname ])
