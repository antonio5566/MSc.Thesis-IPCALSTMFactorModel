clear;
close all;
clc;

%input data
data_name='2020Welch.csv'
data_name_1='fredfactors.mat'
data_name_2='characteristics_data_Sep2023.mat'
Welch=importdata(data_name,',');
DEMEAN=3;
jj=3;
kmax=8;

%% Characteristic Data
load (data_name_2);
T=size(char,1);
charsize=size(char,3);
date_C=yrmo;
char_median= squeeze(nanmedian(char,1));
tcode_char = [5,6,5,5,1,5,5,2,2,2,2,5,5,5,5,5,5,5,5,6,5,5,2,5,5,5,5,5,5,5,5,5,5,5,5,1]
%%
rawdata=Welch.data
% date,index,D12,E12,bm,tbl,AAA,BAA,lty,ntis,Rfree,infl,ltr,corpr,svar,csp
%% Default yield spread (dfy)
dfy=rawdata(:,8)-rawdata(:,7);
%% Term Spread (tms)
tms=rawdata(:,9)-rawdata(:,6);
%% Treasury-bill rate (tbl)
tbl=rawdata(:,6);
%% Book-to-market ratio (bm)
bm=rawdata(:,5);
%% Net equity Expansion (ntis)
ntis=rawdata(:,10);
%% Stock variance (svar) 
svar=rawdata(:,15);
%% Divident-price ratio (dp)
dp= log(rawdata(:,3))-log(rawdata(:,2));
%% Earnings-price ratio (Ep)
ep= log(rawdata(:,4))-log(rawdata(:,2));

%% Combine the data 
date_W=rawdata(:,1);
Welchdata=[dp,ep,bm,ntis,tbl,tms,dfy,svar];

%% Transformation code
tcode_W=[2,2,5,2,2,1,2,1];

%% Announce Time series
load(data_name_1);

%% Combine t_code
t_code=[tcode_F,tcode_char,tcode_W];

%% Match dataset

[loc1 ,inddateex_1] = ismember(date_W,date_F);
[loc2 ,index_2] = ismember(date_W,date_C);
[loc3 ,index_3] = ismember(date_F,date_C);
%% Combine all data
first_occurr_F = find(loc3==1, 1, 'first');
first_occurr_W = find(loc2==1, 1, 'first');
last_occurr_F = find(loc3==1, 1, 'last');
last_occurr_W = find(loc2==1, 1, 'last');

data_F=rawdata_F(first_occurr_F:last_occurr_F,:);
data_W=Welchdata(first_occurr_W:last_occurr_W,:);
data_all=[data_F,char_median,data_W];
T=size(date_C,1);
data_all=data_all(1:T,:);

%% PROCESS DATA
yt=prepare_missing(data_all,t_code);
yt=yt(3:T,:);
date_C=date_C(3:T,:);
[data,n]=remove_outliers(yt);
%% missing value


data_1=fillmissing(data,'nearest');
data-data_1
%% Data name
v_name_W= {'dp','ep','bm','ntis','tbl','tms','dfy','svar'};
v_name=[v_name_F,charnames,v_name_W]
A=[v_name;num2cell(data_1)];

v_name_d={'date'}
B=[v_name_d;num2cell(date_C)];
data_macro=[B,A]


%%
save ("macrodata.mat","data","n","v_name","date_C");
%writecell(data_macro,'macro.csv')


