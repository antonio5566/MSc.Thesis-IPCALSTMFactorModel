clear all
dataname = 'characteristics_data_Nov2023'
dataname_1 = 'macrodata'
load([dataname]);
load([dataname_1]);

%% get median of return to capture the trend
[data_ret,n]=remove_outliers(ret);
medianret= nanmedian(data_ret,1)

%% optimal lag 


[t, variables] = size(data);

maxLag = 12;  

optimalLags = zeros(1, variables);

for variable = 1:variables
    aic_values = zeros(1, maxLag);
    
    for lag = 1:maxLag

        laggedVariable = lagmatrix(data(:, variable), 1:lag);
        % Fit an ARIMA model
        model = arima('ARLags', 1:lag);  
        fit = estimate(model, data(:, variable));
        residuals = infer(fit, data(:, variable));
        % Calculate the log-likelihood
        logLikelihood = -0.5 * (t * log(2 * pi) + sum(residuals.^2));
        numParams = length(fit.AR) + length(fit.MA) + 1;
        aic_values(lag) = 2 * numParams - 2 * logLikelihood;
    end

    [~, optimalLag] = min(aic_values);
    optimalLags(variable) = optimalLag;
end
%%
% ... (your existing code for determining optimal lags)

% Calculate MSE for each variable using optimal lags
mse_values = zeros(1, variables);

for variable = 1:variables
    optimalLag = optimalLags(variable);
    
    % Use the optimal lag to fit the ARIMA model
    laggedVariable = lagmatrix(data(:, variable), 1:optimalLag);
    model = arima('ARLags', 1:optimalLag);
    fit = estimate(model, data(:, variable));
    
    % Infer residuals
    residuals = infer(fit, data(:, variable));
    
    % Calculate MSE
    mse_values(variable) = mean(residuals.^2);
    
    % Display MSE for each variable
    fprintf('MSE for Variable %d: %.4f\n', variable, mse_values(variable));
end

% Display overall MSE
overall_mse = mean(mse_values);
fprintf('Overall MSE: %.4f\n', overall_mse);

v_name_lag = cell(size(v_name));

%%
data_1 = fillmissing(data, 'movmean', 3)
for i = 1:variables
    lag=optimalLags(i);
    laggedmacro(:, i) = [NaN(lag, 1); data(1:end-lag, i)];
    %laggedmacro(:, i) = NaN(size(data, 1), 1);  % Initialize the lagged column
    
    %for j = lag:-1:1
        %laggedmacro(j:end, i) = data(1:end-j+1, i);
        %v_name_lag{end+1} = [v_name{i}, sprintf('_lag%d', j)];
    %end
end
%%
mf_1= [data]
%mf_1= [laggedmacro]
%%
v_name_lag = cell(size(v_name));
%Add "lag_3" after every cell
for i = 1:length(v_name)
    v_name_lag{i} = [v_name{i}, '_lag'];
end

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

