%% Set estimation parameters and data choice
clear all
%mkdir Empirical_Result/macrofactor
%mkdir Empirical_Result/significant
%%

    clearvars -except K
macrochoices = {'macrosLASSOtimestep12'};


dataname = 'IPCADATA_FNW36_RNKDMN_CON'
% Load data
load([dataname]);
Nts = sum(LOC);
% als_opt
als_opt.MaxIterations       = 5000;
als_opt.Tolerance           = 1e-6;


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
%csv_in_3='macro_factor_step12_raw.csv'
%macro_factor_3=importdata(csv_in_3);
%mfactors_3=macro_factor_3.data(1:end,1:4)
%%
%csv_in_4='macro_factor_step18_raw.csv'
%macro_factor_4=importdata(csv_in_4);
%mfactors_4=macro_factor_4.data(1:end,1:4)
%%
%csv_in_5='macro_factor_step24_raw.csv'
%macro_factor_5=importdata(csv_in_5);
%mfactors_5=macro_factor_5.data(1:end,1:4)
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
%%

mf = [mfactor_0(loc1) mfactor_1(loc1) mfactor_2(loc1) mfactor_3(loc1) mfactor_4(loc1) mfactor_5(loc1) mfactor_6(loc1) mfactor_7(loc1) mfactor_8(loc1) mfactor_9(loc1) mfactor_10(loc1) mfactor_11(loc1)]';

%%
[data_ret,n]=remove_outliers(xret);
medianret= nanmedian(data_ret,1)
ret= medianret'


%% LASSO
[Beta, stats] = lasso(mf', ret, 'CV', 5);
lassoPlot(Beta,stats,'PlotType','CV')
%%
Blasso= [Beta(:,stats.Index1SE)];
minMSE= stats.IndexMinMSE;
select = find(Blasso ~= 0);
macro=mf(:,select);
macro=fillmissing(macro,'nearest');
%%
v_name_pre=[v_name,v_name_lag];
v_name=v_name_pre(:,select);
%%
A=[v_name;num2cell(macro)];
%%
v_name_d={'date'}
B=[v_name_d;num2cell(date_C)];
%%
data_macro=[B,A]
writecell(data_macro,'macro_1.csv')
