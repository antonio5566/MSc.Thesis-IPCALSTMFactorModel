clear
close all
clc

csv_in='2015-01.csv';
jj=3;
kmax=8;
dum=importdata(csv_in,',');
series=dum.textdata(1,2:end);
tcode_F=dum.data(1,:);

% Raw data
rawdata_F=dum.data(2:end,:);

% Month/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_month=final_datevec(2);
final_year=final_datevec(1);

% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
t_1=datetime(1959,1,1);
t_2=datetime(final_year,final_month,1)
t1=(t_1:calmonths(1):t_2)'
date_F=datetime(t1,'Format','yyyyMM')
% T = number of months in sample
T=size(date_F,1);
rawdata_F=rawdata_F(1:T,:);

date_f_1=string(date_F,'yyyyMM')
date_F=double(date_f_1)
save ("fredfactors.mat","rawdata_F","date_F","tcode_F")

%%
% =========================================================================
% PART 2: PROCESS DATA

% Transform raw data to be stationary using auxiliary function
% prepare_missing()
yt=prepare_missing(rawdata_F,tcode_F);

% Reduce sample to usable dates: remove first two months because some
% series have been first differenced
yt=yt(3:T,:);
date_F=date_F(3:T,:);

% Remove outliers using auxiliary function remove_outliers(); see function
% or readme.txt for definition of outliers
%   data = matrix of transformed series with outliers removed
%   n = number of outliers removed from each series
[data_F,n]=remove_outliers(yt);
v_name_F=dum.textdata(1,2:end)
save ("fredfactors.mat","rawdata_F","date_F","v_name_F","tcode_F")
