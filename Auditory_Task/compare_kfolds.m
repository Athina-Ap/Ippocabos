% This function computes the mean and standard deviation of the R-squared
% computed for each of the folds in the k-fold cross validation for the
% PGAM. The aim is to understand how much the results from each of the
% folds vary. 

function [r2_full_mean, r2_full_std, r2_reduced_mean, r2_reduced_std] ...
    = compare_kfolds(varargin)

p = inputParser;
addParameter(p,'basepath','F:\PGAM\fits_new\k_fold_tests',@isstr);
parse(p,varargin{:});
basepath = p.Results.basepath;

% list all reduced and full .csv files
% FIXME: the neurons are not in order. 
filelist_full = dir(fullfile(basepath, '*full.csv'));
filelist_reduced = dir(fullfile(basepath, '*reduced.csv'));

neurons_full = zeros(length(filelist_full(:,1)), 1);
neurons_reduced = zeros(length(filelist_reduced(:,1)), 1);

filename_full = {filelist_full.name};
filename_reduced = {filelist_reduced.name};

for i = 1:length(filename_full(:))
    neurons_full(i,1) = str2double(regexp(filename_full{i}, '\d+', 'match'));
end
   
for i = 1:length(filename_reduced(:))
    neurons_reduced(i,1) = str2double(regexp(filename_reduced{i}, '\d+', 'match'));
end

neurons = intersect(neurons_full, neurons_reduced);

% Array of r2 values for each fold 
r2_full = zeros(5, length(neurons));
r2_full_mean = zeros(1, length(neurons));
r2_full_std = zeros(1, length(neurons));

r2_reduced = zeros(5, length(neurons));
r2_reduced_mean = zeros(1, length(neurons));
r2_reduced_std = zeros(1, length(neurons));

for i = 1:length(neurons)   
    table = readtable(fullfile(basepath, filename_full{neurons(i)}), ...
        'readvariablenames', false);

    r2_full(:,i) = table.Var2(2:end,1);
    r2_full_mean(i) = mean(r2_full(:,i));
    r2_full_std(i) = std(r2_full(:,i));
end

for i = 1:length(neurons)
    table = readtable(fullfile(basepath, filename_reduced{i}), ...
        'readvariablenames', false);

    r2_reduced(:,i) = table.Var2(2:end,1);
    r2_reduced_mean(i) = mean(r2_reduced(:,i));
    r2_reduced_std(i) = std(r2_reduced(:,i));
end

figure;
hold on;
colors = ['b', 'r'];

p1 = scatter(neurons(:)', r2_full(:,:), 8, 'bo', 'MarkerFaceColor', ...
        colors(1), 'DisplayName', 'Full');
p2 = scatter(neurons(:)', r2_reduced(:,:), 8, 'ro', 'MarkerFaceColor', ...
        colors(2), 'DisplayName', 'Reduced');

ylim([0,0.6])
legend([p1(1) p2(1)])
xlabel('Neuron');
ylabel('Pseudo R-squared');
title(['Pseudo R-squared values for the different folds for the PGAM ' ...
    'cross-validation for a sample of neurons.'])


end