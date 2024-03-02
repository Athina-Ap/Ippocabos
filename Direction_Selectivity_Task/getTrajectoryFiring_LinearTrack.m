function getTrajectoryFiring_LinearTrack(placeCells, varargin)

% This is an adaptation of the function plotting7port.m (Ipshita Zutshi).

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveLoc',[],@isstr);
addParameter(p,'plotFig','true',@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
saveLoc = p.Results.saveLoc;
plotFig = p.Results.plotFig;

if isempty(saveLoc)
    saveLoc = strcat(basepath,'\TrajectoryMaps');
    if ~isfolder('TrajectoryMaps')
        mkdir('TrajectoryMaps')
    end    
end
%% Deal with inputs
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikeData.cellinfo.mat']))
    file = dir([basepath filesep '*.spikeData.cellinfo.mat']);
    load(file.name);
end

if ~isempty(dir([basepath filesep '*.rateMapsAvg.cellinfo.mat']))
    file = dir([basepath filesep '*.rateMapsAvg.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

if strcmp(behavTrials.start, 'left')
    labels{1} = 'right';
    labels{2} = 'left';
else
    labels{1} = 'left';
    labels{2} = 'right';
end

if plotFig 
    for cellNum = 1:length(placeCells)
        for pf = 1:(size(behavTrials.timestamps,1)-1)    
            [idx{1}] = InIntervals(tracking.timestamps, behavTrials.timestamps(pf,:));
            positions.(labels{1}){pf} = [tracking.position.x(idx{1}) tracking.position.y(idx{1})];   
            [idx{2}] = InIntervals(tracking.timestamps, [behavTrials.timestamps(pf,2) behavTrials.timestamps(pf+1,1)]);
            positions.(labels{2}){pf} = [tracking.position.x(idx{2}) tracking.position.y(idx{2})];   
        end
          
        figure
        set(gcf,'Color','w')
        set(gcf,'Position',[2050 181 1585 762])
        plot(tracking.timestamps,tracking.position.y,'Color',[160/243 160/243 160/243], 'LineWidth',1.2)
        hold on
        scatter(tracking.timestamps(spikeData.posIdx{placeCells(cellNum)}),tracking.position.y(spikeData.posIdx{placeCells(cellNum)}),8,'r','filled')
        ylim([0 500])
        xlabel('Time(s)')
        ylabel('Position on track (px)')
        saveas(gcf,[saveLoc,filesep ,'cell_', num2str(placeCells(cellNum)),'.png'],'png');
        saveas(gcf,[saveLoc,filesep ,'cell_', num2str(placeCells(cellNum)),'.fig'],'fig');
        close all 
    end
end
end

