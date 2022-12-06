function [new_metrics_table, varNames] = Metrics_Compute(r,x,Fs,MetricsToCompute)
% INPUT PARAMETERS
%   x                 -Processed or enhanced signal
%   r                 -Reference signal used to compute the metrics
%   Fs                -Samplig frequency
%   MetricsToCompute  -cell array defining the metrics to compute 'stoi',
%                       'sd'
% OUTPUT PARAMETERS
%   new_metrics_table - Table with metrics
%   varNames          - cell array with Metrics names


metrics_value = {};
j=1;
varNames ={};
for i =1:length(MetricsToCompute)
    switch MetricsToCompute{i}  
        case {'stoi','STOI'}
            stoi_est = stoi(r,x, Fs);
            metrics_value{j} = stoi_est;
            varNames{j} = 'STOI';
        case {'sd','SD'}
            SD_scores  = SD_measure_fbcancellation(r,x, Fs, length(r)/Fs); % Distortion measure
            metrics_value{j} = SD_scores(1);
            metrics_value{j+1} = SD_scores(2);
            varNames{j} = 'mean SD';
            varNames{j+1} = 'max SD';
            j = j+1;
    end
    j=j+1;
end

metrics_table = cell2table(metrics_value,'VariableNames',varNames);
n_decimal=2;
new_metrics_table = varfun(@(x) num2str(x, ['%' sprintf('.%df', n_decimal)]), metrics_table);
% preserve the variable names and the row names in the original table
new_metrics_table.Properties.VariableNames = metrics_table.Properties.VariableNames;
