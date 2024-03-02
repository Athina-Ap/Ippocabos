function [directionStats] = getPlaceFieldDrift(varargin)
% Drift/stability analysis for place fields.
% This is done by splitting the session into the first and the second half
% and correlating the two halves. 

close all 

p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;

%% Load data
% Get session info
basename = bz_BasenameFromBasepath(basepath);
sessionInfo = load([basepath filesep basename '.sessionInfo.mat']);

% Get place fields
if ~isempty([basepath filesep basename 'placeFields.cellinfo.mat'])
    load([basepath filesep basename '.placeFields.cellinfo.mat']);
else
    disp('No place fields in this directory.')
    exit
end

% Get behavior
if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

% Get merge points
if ~isempty(dir([basepath filesep '*.MergePoints.events.mat'])) 
    file = dir([basepath filesep '*.MergePoints.events.mat']);
    load(file(1).name);
end

% Get firing maps 
if ~isempty([basepath filesep basename 'firingMapsTrial.cellinfo.mat'])
    load([basepath filesep basename '.firingMapsTrial.cellinfo.mat']);
else
    disp('No rate maps in the this directory.')
    exit 
end

if exist('MergePoints') & length(MergePoints.foldernames) > 1
    if ~isempty([basepath filesep basename 'firingMapsCond.cellinfo.mat'])
        load([basepath filesep basename '.firingMapsCond.cellinfo.mat']);
    end
end

% Get direction stats
if ~isempty(dir([basepath filesep '*.directionStats.cellinfo.mat'])) 
    file = dir([basepath filesep '*.directionStats.cellinfo.mat']);
    load(file(1).name);
end

%% Initialize 
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    labels{1} = 'right';
    labels{2} = 'left';

    conditions{1} = '1port';
    conditions{2} = '2ports';
else
    if strcmp(placeFieldStats.params.start, 'left')
        labels{1} = 'right';
        labels{2} = 'left';
    else
        labels{1} = 'left';
        labels{2} = 'right';
    end
end 

%% Calculate drift/stability
numtrials = length(firingMaps.right.rateMaps{1,1});

if exist('MergePoints') & length(MergePoints.foldernames) > 1
    
    % Split trials in half
    mergeidx = find(behavTrials.timestamps(:,1) > MergePoints.timestamps(1,2), 1,'first');

    sessLen{1} = 1:mergeidx;
    sessLen{2} = mergeidx:numtrials;
    
    for cond = 1:length(conditions)
        placeCells{cond} = unique([directionStats{cond}.placeCells.(labels{1}), ...
            directionStats{cond}.placeCells.(labels{2})]);
        first_half{cond} = 1:ceil(length(sessLen{cond})/2);
        second_half{cond} = ceil(length(sessLen{cond})/2)+1:length(sessLen{cond});
    end

    % Group firing data for first and second halves and find their means
    for cond = 1:length(conditions)
        for pf = 1:length(placeCells{cond})
            for ll = 1:length(labels)
                datamat.(labels{ll}){cond}{1} = []; % first half
                datamat.(labels{ll}){cond}{2} = []; % second half

                for kk = first_half{cond}
                    datamat.(labels{ll}){cond}{1} = [datamat.(labels{ll}){cond}{1}; ...
                        firingMaps.(labels{ll}).rateMaps{placeCells{cond}(pf)}{kk}];
                end
                datamat_mean.(labels{ll}){cond}{1}{pf} = mean(datamat.(labels{ll}){cond}{1}, 1);

                for kk = second_half{cond}
                    datamat.(labels{ll}){cond}{2} = [datamat.(labels{ll}){cond}{2}; ...
                        firingMaps.(labels{ll}).rateMaps{placeCells{cond}(pf)}{kk}];
                end
                datamat_mean.(labels{ll}){cond}{2}{pf} = mean(datamat.(labels{ll}){cond}{2}, 1);
            end
        end
    end

    % Correlate first and second halves 
    for cond = 1:length(conditions)
        for ll = 1:length(labels)
            corrMap{cond}{ll} = zeros(length(placeCells{cond}),1);
    
            for i = 1:length(placeCells{cond})
                corr = corrcoef(datamat_mean.(labels{ll}){cond}{1}{i}', ...
                    datamat_mean.(labels{ll}){cond}{2}{i}','rows','complete');
                corrMap{cond}{ll}(i) = corr(1,2);
                directionStats{cond}.corr{ll}(i) = corrMap{cond}{ll}(i);
            end
            corrMap_cellsMean{cond}{ll} = nanmean(corrMap{cond}{ll}, 1);
            directionStats{cond}.corrCellsMean{ll} = corrMap_cellsMean{cond}{ll};
        end
        % Final correlation is the mean between two running directions
        corrMap_dirsMean{cond} = mean(cell2mat(corrMap_cellsMean{cond}));
        directionStats{cond}.corrDirsMean = corrMap_dirsMean{cond};
    end

else
    placeCells = unique([directionStats.placeCells.(labels{1}), ...
        directionStats.placeCells.(labels{2})]);

    % Split trials in half
    first_half = 1:ceil(numtrials/2);
    second_half = ceil(numtrials/2)+1:numtrials;
    
    % Group firing data for first and second halves and find their means
    for pf = 1:length(placeCells)  
        for ll = 1:length(labels)
            datamat.(labels{ll}){1} = []; % first half
            datamat.(labels{ll}){2} = []; % second half

            for kk = first_half
                datamat.(labels{ll}){1} = [datamat.(labels{ll}){1}; ...
                    firingMaps.(labels{ll}).rateMaps{placeCells(pf)}{kk}];
            end
            datamat_mean.(labels{ll}){1}{pf} = mean(datamat.(labels{ll}){1}, 1);

            for kk = second_half
                datamat.(labels{ll}){2} = [datamat.(labels{ll}){2}; ...
                    firingMaps.(labels{ll}).rateMaps{placeCells(pf)}{kk}];
            end
            datamat_mean.(labels{ll}){2}{pf} = mean(datamat.(labels{ll}){2}, 1);
        end
        
    end

    % Correlate first and second halves 
    for ll = 1:length(labels)
        corrMap{ll} = zeros(length(placeCells),1);

        for i = 1:length(placeCells)
            corr = corrcoef(datamat_mean.(labels{ll}){1}{i}', ...
                datamat_mean.(labels{ll}){2}{i}','rows','complete');
            corrMap{ll}(i) = corr(1,2);
            directionStats.corr{ll}(i) = corrMap{ll}(i);
        end
        corrMap_cellsMean{ll} = nanmean(corrMap{ll}, 1);
        directionStats.corrCellsMean{ll} = corrMap_cellsMean{ll};
    end
    % Final correlation is the mean between two running directions
    corrMap_dirsMean = mean(cell2mat(corrMap_cellsMean));
    directionStats.corrDirsMean = corrMap_dirsMean;
    
end
                
%% Save output
if saveMat
   save([basepath,filesep,basename '.directionStats.cellinfo.mat'],'directionStats'); 
end

end