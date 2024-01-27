clear
clc
data = importdata('characteristics_data_feb2017.csv');
%%
mkdir Empirical_Result
% The first 4 columns of data.textdata are not data. Col 2 is the year, col 3 the year

yr = cellfun(@str2num,data.textdata(2:end,2),'un',true);
mo = cellfun(@str2num,data.textdata(2:end,3),'un',true);
bigyrmo = yr*100+mo;

bigpermno = data.data(:,1);

yrmo = unique(bigyrmo);
permno = unique(bigpermno);

T = length(yrmo);
N = length(permno);

textdatacol = [];
for j=1:size(data.textdata,2);textdatacol{j} = regexprep(data.textdata{1,j},'"','');end

%%
neuhierlchoices = {...
    'lme' 'lturnover' 'spread_mean' ... mismatch with notation spread
    'cum_return_1_0' 'cum_return_12_2' 'cum_return_12_7' ...
    'cum_return_36_13' 'beta' 'idio_vol' 'beme' ...
    'at' 'ato' 'c' 'cto' 'd2a' 'dpi2a' 'e2p' 'fc2y' 'free_cf' 'investment' 'lev' 'noa' 'oa' 'ol' 'pcm' ...
    'pm' 'prof' 'q' 'rel_to_high_price' 'rna' 'roa' 'roe' 's2p' 'sga2m' 'suv' ...
    'a2me'};

[~,~,charloc] = intersect(neuhierlchoices,textdatacol);
[~,~,retloc]  = intersect('ret',textdatacol);

char = nan(N,T,length(charloc));
charnames = textdatacol(charloc);
ret   = nan(N,T);
%%
for n=1:N
    nloc = find(bigpermno==permno(n));
    tloc = find(ismember(yrmo,bigyrmo(nloc)));
    ret(n,tloc) = data.data(nloc,retloc-4);
    char(n,tloc,:) = data.data(nloc,charloc-4);
end

%%
save('characteristics_data_Nov2023.mat','ret','char','charnames','yrmo','permno','-v7.3')