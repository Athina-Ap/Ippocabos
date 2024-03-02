% This function calculates the R-squared for the PGAM fits for all neurons
% in a session. The aim is to compare the R-squared with different models
% e.g. when the knots for a variable are 10 vs. when they are 20. 

function [r2] = compare_r2(varargin)

p = inputParser;
addParameter(p,'basepath','F:\PGAM\fits_new',@isstr);
parse(p,varargin{:});
basepath = p.Results.basepath;

% list all .csv files
% FIXME: the neurons are not in order. 
filelist = dir(fullfile(basepath, '*fit.csv'));
filelist_highDens = dir(fullfile(basepath, '*fit_highDens.csv'));

% Array of r2 values for each cell
r2 = zeros(1, length(filelist(:,1)));
r2_highDens = zeros(1, length(filelist_highDens(:,1)));

for i = 1:length(filelist(:,1))
    filename = filelist(i,1).name;

    table = readtable(fullfile(basepath, filename), 'readvariablenames', false);

    % r2 is the same for all variables - only select on row to read from
    r2(i) = table.Var6(2);
end

for i = 1:length(filelist_highDens(:,1))
    filename_highDens = filelist_highDens(i,1).name;

    table_highDens = readtable(fullfile(basepath, filename_highDens), 'readvariablenames', false);

    % r2 is the same for all variables - only select on row to read from
    r2_highDens(i) = table_highDens.Var6(2);
end

% In case the two have different sizes (temporary).
% Create scatter plot of r2. 
if length(r2) > length(r2_highDens)
    r2_subset = r2(1:length(r2_highDens));

    figure;
    scatter(r2_subset, r2_highDens, 'filled');
    xlim([0,1])
    ylim([0,1])
    hold on 
    plot([0,1], [0,1]);
    hold off

elseif length(r2) < length(r2_highDens)
    r2_highDens_subset = r2_highDens(1:length(r2));

    figure;
    scatter(r2, r2_highDens_subset, 'filled');
    xlim([0,1])
    ylim([0,1])
    hold on
    plot([0,1], [0,1]);
    hold off

else 
    r2_subset = r2;

    figure;
    scatter(r2_subset, r2_highDens, 'filled');
    xlim([0,1])
    ylim([0,1])
    hold on 
    plot([0,1], [0,1]);
    hold off
end

end