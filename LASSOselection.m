clear all
dataname = 'characteristics_data_Nov2023'
dataname_1 = 'macrodata'
load([dataname]);
load([dataname_1]);

%% get median of return to capture the trend
[data_ret,n]=remove_outliers(ret);
medianret= nanmedian(data_ret,1)

%% standardization of macro
mf_2 = tiedrank(mf_1);
mf = (mf_2-1)./(max(mf_2)-1);
mf = mf-0.5;

%%
[date,loc1,loc2] = intersect(yrmo,date_C);
ret = medianret (:,loc1)'

%% LASSO
[Beta, stats] = lasso(mf, ret, 'CV', 5);
lassoPlot(Beta,stats,'PlotType','CV')
%%
Blasso= [Beta(:,stats.IndexMinMSE)];
minMSE= stats.IndexMinMSE;
select = find(Blasso ~= 0);
macro=mf_1(:,select);
macro=fillmissing(macro,'nearest');
%%
v_name_pre=[v_name];
v_name_2 =v_name_pre(:,select);

%%
A=[v_name_2;num2cell(macro)];
%%
v_name_d={'date'}
B=[v_name_d;num2cell(date_C)];
%%
data_macro=[B,A]
writecell(data_macro,'macro_1.csv')
%%

