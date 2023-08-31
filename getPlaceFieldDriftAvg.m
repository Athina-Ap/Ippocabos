function [outputArg1,outputArg2] = getPlaceFieldDriftAvg(varargin)
% Drift/stability analysis for place fields.
% This is done by splitting the session into the first and the second half
% and correlating the two halves. 

close all 

p=inputParser;
addParameter(p,'root',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'analyze2sess',true,@islogical);

parse(p,varargin{:});
root = p.Results.root;
saveMat = p.Results.saveMat;
analyze2sess = p.Results.analyze2sess;

%% Load data and pre-define variables
corr_dirsMean = cell(3, 3); % conditions, max sessions/condition

% Get average direction stats from all sessions
if ~isempty(dir([root filesep 'directionStatsAvg.cellinfo.mat'])) 
    file = dir([root filesep 'directionStatsAvg.cellinfo.mat']);
    load(file(1).name);
end

for c = 1:3
    corr_sessMean{c} = [];

    if c == 1
        sessions = ["IZ50/IZ50_230621_sess9", "IZ51/IZ51_230621_sess9", ...
            "IZ50/IZ50_230623_sess12"];
    elseif c == 2
        sessions = ["IZ50/IZ50_230621_sess10", "IZ51/IZ51_230621_sess10", ...
            "IZ50/IZ50_230623_sess12"];
    else
        sessions = ["IZ50/IZ50_230622_sess11", "IZ51/IZ51_230622_sess11"];
    end

    for s = 1:length(sessions)
        cd(sessions(s))
        basepath = pwd;

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
        
        % Calculate mean of correlation for each cell across directions   
        if exist("MergePoints") & length(MergePoints.foldernames) > 1
            if c == 1
                for i = 1:length(directionStats{2}.corr{1})
                    corr_dirsMean{c}{s}(i) = mean([directionStats{2}.corr{1}(i), ...
                        directionStats{2}.corr{2}(i)], 2);
                end
            elseif c == 2
                for i = 1:length(directionStats{1}.corr{1})
                    corr_dirsMean{c}{s}(i) = mean([directionStats{1}.corr{1}(i), ...
                        directionStats{1}.corr{2}(i)], 2);
                end
            end
            
        else
            for i = 1:length(directionStats.corr{1})
                corr_dirsMean{c}{s}(i) = mean([directionStats.corr{1}(i), ...
                        directionStats.corr{2}(i)], 2);
            end
        end

        % Concatenate average correlations across session
        corr_sessMean{c} = [corr_sessMean{c}; corr_dirsMean{c}{s}'];
        cd(root)
    end 
end

%% Statistical analysis 
colors = [[0, 0.2471, 0.3608]; ...
        [0.5843, 0.3176, 0.5882]; ...
        [1, 0.6510, 0]];
data = corr_sessMean;
% group = [ones(length(corr_sessMean{1}),1); ones(length(corr_sessMean{2}),1)*2; ...
%         ones(length(corr_sessMean{1}),1)*3]; 
    
directionStatsAvg.stability.stats_all = groupStats(data, [], 'color',colors, ...
    'plotType','boxplot','labelSummary',false,'sigStarTest','KW');

set(gcf,'Position',[300 100 700 600]);
ylabel('Corr. Coefficient', FontSize=16)
xticks([1,2,3])
xticklabels(["2 ports", "1 port", "walls off"])
xticklabel = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xticklabel, 'FontSize', 16)

saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability.fig'],'fig');
close(figure(1))

% Separate out by condition 
colors1 = [[0, 0.2471, 0.3608]; ...
        [1, 0.6510, 0]];
colors2 = [[0, 0.2471, 0.3608]; ...
        [0.5843, 0.3176, 0.5882]];
data1 = {corr_sessMean{1}; corr_sessMean{3}};
data2 = {corr_sessMean{1}; corr_sessMean{2}};

directionStatsAvg.stability.stats_2portswallsoff = groupStats(data1, [], 'color',colors1, ...
    'plotType','boxplot','labelSummary',false,'sigStarTest','KW');

set(gcf,'Position',[300 100 700 600]);
ylabel('Corr. Coefficient', FontSize=16)
xticks([1,2])
ylim([-1,1.2])
yticks([-1,0,1])
xticklabels(["2 ports", "walls off"])
xticklabel = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xticklabel, 'FontSize', 16)
saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability_2portsWallsOff.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability_2portsWallsOff.fig'],'fig');

directionStatsAvg.stability.stats_2ports1port = groupStats(data2, [], 'color',colors2, ...
    'plotType','boxplot','labelSummary',false,'sigStarTest','KW');

set(gcf,'Position',[300 100 700 600]);
ylabel('Corr. Coefficient', FontSize=16)
xticks([1,2,3])
ylim([-1,1.2])
yticks([-1,0,1])
xticklabels(["2 ports", "1 port"])
xticklabel = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xticklabel, 'FontSize', 16)
saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability_2ports1port.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability_2ports1port.fig'],'fig');

%% Make extra plots for double sessions
if analyze2sess
    colors = [[0.5843, 0.3176, 0.5882]; [0, 0.2471, 0.3608]];
    data2sess = {corr_dirsMean{1}{3}, corr_dirsMean{2}{3}};

    directionStatsAvg.stability2sess.stats = groupStats(data2sess, [], 'color',colors, ...
        'plotType','boxplot','labelSummary',false,'sigStarTest','KW');

    set(gcf,'Position',[300 100 700 600]);
    ylabel('Corr. Coefficient', FontSize=16)
    xticks([1,2,3])
    xticklabels(["1 port", "2 ports"])
    xticklabel = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', xticklabel, 'FontSize', 16)

    saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability_2sess.png'],'png');
    saveas(gcf,[root filesep 'FiringMapAvg\correlations',filesep,'stability_2sess.fig'],'fig');
    close(figure(1))
end
    
%% Save output
if saveMat
   save([root,filesep 'directionStatsAvg.cellinfo.mat'],'directionStatsAvg'); 
end
close all
end