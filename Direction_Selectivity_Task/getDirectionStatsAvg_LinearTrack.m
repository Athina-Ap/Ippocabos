function [directionStatsAvg] = getDirectionStatsAvg_LinearTrack(varargin)

% (1) Plot rate maps of place cells from all sessions.
% (2) Perform correlation analyses with all sessions.

close all 

p=inputParser;
addParameter(p,'root',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'thresFR',0.5,@isnumeric);

parse(p,varargin{:});
root = p.Results.root;
saveMat = p.Results.saveMat;
thresFR = p.Results.thresFR;

%% Load data and pre-define variables
datamat_cells{1} = cell(3,1);
datamat_cells{2} = cell(3,1);

direction{1} = 'right';
direction{2} = 'left';

% Get average direction stats from all sessions
if ~isempty(dir([root filesep 'directionStatsAvg.cellinfo.mat'])) 
    file = dir([root filesep 'directionStatsAvg.cellinfo.mat']);
    load(file(1).name);
end

% c = 1 [2 ports], c = 2 [1 port], c = 3 [walls off]
for c = 1:3

    if c == 1
        sessions = ["IZ50/IZ50_230621_sess9", "IZ51/IZ51_230621_sess9", ...
            "IZ50/IZ50_230623_sess12"];
    elseif c == 2
        sessions = ["IZ50/IZ50_230621_sess10", "IZ51/IZ51_230621_sess10", ...
            "IZ50/IZ50_230623_sess12"];
    else
        sessions = ["IZ50/IZ50_230622_sess11", "IZ51/IZ51_230622_sess11"];
    end
        
    if ~isfolder([root filesep 'FiringMapAvg/PlaceCells'])
        mkdir([root filesep 'FiringMapAvg/PlaceCells']);
    end
    if ~isfolder([root filesep 'FiringMapAvg/correlations'])
        mkdir([root filesep 'FiringMapAvg/correlations']);
    end

    for s = 1:length(sessions)
        cd(sessions(s))
        basepath = pwd;

        %% Load data
        % Get session info
        basename = bz_BasenameFromBasepath(basepath);
        load([basepath filesep basename '.sessionInfo.mat']);
        
        % Get firing maps 
        if ~isempty([basepath filesep basename '.firingMapsTrial.cellinfo.mat'])
            load([basepath filesep basename '.firingMapsTrial.cellinfo.mat']);
        else
            disp('No rate maps in the this directory.')
            exit 
        end
        if exist([basepath filesep basename '.firingMapsCond.cellinfo.mat'])
            load([basepath filesep basename '.firingMapsCond.cellinfo.mat']);
        end
        
        % Get place fields
        if ~isempty([basepath filesep basename '.placeFields.cellinfo.mat'])
            load([basepath filesep basename '.placeFields.cellinfo.mat']);
        else
            disp('No place fields in this directory.')
            exit
        end
        
        % Get tracking
        if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
            file = dir([basepath filesep '*.Tracking.Behavior.mat']);
            load(file(1).name);
        end
        
        % Get behavior
        if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
            file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
            load(file(1).name);
        end
        
        % Get spikes 
        if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
            file = dir([basepath filesep '*.spikes.cellinfo.mat']);
            load(file.name);
        end

        % Get MergePoints
        if ~isempty(dir([basepath filesep '*.MergePoints.events.mat']))
            file = dir([basepath filesep '*.MergePoints.events.mat']);
            load(file.name);
        end
                
        if exist('MergePoints') & length(MergePoints.foldernames) > 1
            labels{1} = 'right';
            labels{2} = 'left';
        else
            if strcmp(placeFieldStats.params.start, 'left')
                labels{1} = 'right';
                labels{2} = 'left';
            else
                labels{1} = 'left';
                labels{2} = 'right';
            end
        end
    
        %% Get cells with place field(s)
        for d = 1:length(labels)
            placeCells{d} = [];
            if exist('MergePoints') & length(MergePoints.foldernames) > 1
                if c == 1
                    for i = 1:length(placeFieldStats{2}.mapStats)
                        if ~isnan(placeFieldStats{2}.mapStats{i}{d}.x)
                            placeCells{d} = [placeCells{d}, i];
                        end
                    end
                elseif c == 2
                    for i = 1:length(placeFieldStats{1}.mapStats)
                        if ~isnan(placeFieldStats{1}.mapStats{i}{d}.x)
                            placeCells{d} = [placeCells{d}, i];
                        end
                    end
                end
            else
                for i = 1:length(placeFieldStats.mapStats) 
                    if ~isnan(placeFieldStats.mapStats{i}{d}.x)
                        placeCells{d} = [placeCells{d}, i];
                    end
                end   
            end
        end
    
        %% Concatenate data across sessions
     
        % Apply firing rate threshold
%         dt = mean(diff(tracking.timestamps));
%         win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
%         spkData = bz_SpktToSpkmat(spikes,'dt',dt,'win',win);
%         
%         spkMat = spkData.data';
%         timestamps = spkData.timestamps';
%         
%         fr = zeros(1, length(spkMat(:,1)));
%         for i = 1:length(spkMat(:,1))
%             fr(i) = mean(spkMat(i,:)) * 1/dt;
%         end
%         LowFiringCells = find(fr < thresFR);
        
        % Concatenate data of the same condition
        placeCellsAll = unique([placeCells{1}, placeCells{2}]);

        datamat_mean{1} = cell(length(placeCellsAll),1);
        datamat_mean{2} = cell(length(placeCellsAll),1);
    
        for pf = 1:length(placeCellsAll)   
%             if ismember(placeCellsAll(pf), LowFiringCells)
%                 continue
%             end
            datamat{1} = [];
            datamat{2} = [];
            for ll = 1:length(labels)
               if exist('MergePoints') & length(MergePoints.foldernames) > 1
                   mergeidx = find(behavTrials.timestamps(:,1) > MergePoints.timestamps(1,2), 1,'first');
                   if c == 1
                       for kk = mergeidx:length(firingMaps.right.rateMaps{1})
                           datamat{ll} = [datamat{ll}; firingMaps.(labels{ll}).rateMaps{placeCellsAll(pf)}{kk}];
                       end
                   elseif c ==2
                       for kk = 1:mergeidx
                           datamat{ll} = [datamat{ll}; firingMaps.(labels{ll}).rateMaps{placeCellsAll(pf)}{kk}];
                       end
                   end
               else
                   for kk = 1:length(firingMaps.right.rateMaps{1})
                       datamat{ll} = [datamat{ll}; firingMaps.(labels{ll}).rateMaps{placeCellsAll(pf)}{kk}];
                   end
               end
               datamat_mean{ll}{pf} = mean(datamat{ll}, 1); 
            end
        end
    
        for pf = 1:length(placeCellsAll)
            if strcmp(labels{1}, direction{1})
                for ll = 1:length(labels)
                    datamat_cells{ll}{c} = [datamat_cells{ll}{c}; datamat_mean{ll}{pf}];
                end
            else
                for ll = 1:length(labels)
                    od = 3-ll;
                    datamat_cells{ll}{c} = [datamat_cells{ll}{c}; datamat_mean{od}{pf}];
                end
            end
        end
        cd(root)
    end

    %% Normalize firing rate for each cell
    datamat_cells_normAll{c} = zeros(size(datamat_cells{1}{c},1), 2*size(datamat_cells{1}{c},2));
    for pf = 1:size(datamat_cells{1}{c},1)
        datamat_cells_normAll{c}(pf,:) = normalize([datamat_cells{1}{c}(pf,:), datamat_cells{2}{c}(pf,:)], 'range');
        norm_data = datamat_cells_normAll{c}(pf,:);
        datamat_cells_norm{1}{c}(pf,:) = norm_data(1:length(datamat_cells{1}{c}(pf,:)));
        datamat_cells_norm{2}{c}(pf,:) = norm_data(length(datamat_cells{1}{c}(pf,:))+1:end);
    end
end 

%% Plot place maps for the cells with fields only

% Raw rate
figure;
set(gcf,'Renderer','painters')
set(gcf, 'Position', [1922,204,1867,792]);
k = 1;
for c = 1:3
    for ll = 1:length(labels)
        [~,idx] = max(datamat_cells{1}{c}, [], 2);
        [~,sortidx] = sort(idx);
        datamat_cellsSorted{ll} = datamat_cells{ll}{c}(sortidx, :);
    
        subplot(2,6,k)
        imagesc(datamat_cellsSorted{ll});
        cbar = colorbar;
        cbar.Ticks = [ceil(min(datamat_cellsSorted{ll}(:))), round(max(datamat_cellsSorted{ll}(:)))];
        cbar.FontSize = 10;
        xticks([1,size(datamat_cellsSorted{ll},2)])
        xticklabels([0,120])
        yticks([1,size(datamat_cellsSorted{ll},1)])
        yticklabels([1,size(datamat_cellsSorted{ll},1)])
        set(gca, 'FontSize', 14)
        xlabel('Position (cm)')
        ylabel('Neuron')
        k = k + 6;
    end
    k = k - 11;
    for ll = 1:length(labels)
        [~,idx] = max(datamat_cells{2}{c}, [], 2);
        [~,sortidx] = sort(idx);
        datamat_cellsSorted{ll} = datamat_cells{ll}{c}(sortidx, :);
    
        subplot(2,6,k)
        imagesc(datamat_cellsSorted{ll});
        cbar = colorbar;
        cbar = colorbar;
        cbar.Ticks = [ceil(min(datamat_cellsSorted{ll}(:))), round(max(datamat_cellsSorted{ll}(:)))];
        cbar.FontSize = 10;
        xticks([1,size(datamat_cellsSorted{ll},2)])
        xticklabels([0,120])
        yticks([1,size(datamat_cellsSorted{ll},1)])
        yticklabels([1,size(datamat_cellsSorted{ll},1)])
        set(gca, 'FontSize', 14)
        xlabel('Position (cm)')
        ylabel('Neuron')
        k = k + 6;
    end
    k = k - 11;
end
annotation(figure(1),'textbox',...
    [0.0583824317086234 0.732585859488961 0.0372254944974622 0.0429292920261922],...
    'String',{'right'},'LineStyle','none','FontWeight','bold','FontSize',14);
annotation(figure(1),'textbox',...
    [0.0589180503481521 0.255313132216234 0.0297268337675385 0.0429292920261923],...
    'String',{'left'},'LineStyle','none','FontWeight','bold','FontSize',14);
annotation(figure(1),'textbox',...
    [0.215318693090519 0.943444445347551 0.0484734855923477 0.0429292920261922],...
    'String',{'2 ports'},'LineStyle','none', 'FontWeight','bold','FontSize',14, ...
    'Color',[0, 0.2471, 0.3608]);
annotation(figure(1),'textbox',...
    [0.487412961971076 0.938393940297045 0.0425816807331219 0.0429292920261922],...
    'String',{'1 port'},'LineStyle','none','FontWeight','bold','FontSize',14, ...
    'Color',[0.5843, 0.3176, 0.5882]);
annotation(figure(1),'textbox',...
    [0.750937332619175 0.942181819084922 0.0549009090751394 0.0429292920261922],...
    'String',{'walls off'},'LineStyle','none','FontWeight','bold','FontSize',14, ...
    'Color',[1, 0.6510, 0]);

saveas(gcf,[root filesep 'FiringMapAvg/PlaceCells',filesep,'heatmaps_raw.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/PlaceCells',filesep,'heatmaps_raw.fig'],'fig');
close(figure(1))

% Normalized rate
figure;
set(gcf,'Renderer','painters')
set(gcf, 'Position', [1922,204,1867,792]);
k = 1;
for c = 1:3
    for ll = 1:length(labels)
        [~,idx] = max(datamat_cells_norm{1}{c}, [], 2);
        [~,sortidx] = sort(idx);
        datamat_cellsSorted_norm{ll} = datamat_cells_norm{ll}{c}(sortidx, :);
    
        subplot(2,6,k)
        
        imagesc(datamat_cellsSorted_norm{ll});
        xticks([1,size(datamat_cellsSorted_norm{ll},2)])
        xticklabels([0,120])
        yticks([1,size(datamat_cellsSorted_norm{ll},1)])
        yticklabels([1,size(datamat_cellsSorted_norm{ll},1)])
        set(gca, 'FontSize', 14)
        xlabel('Position (cm)')
        if k == 1 || k == 7
            disp('true')
            ylabel('Neuron')
        end
        k = k + 6;
    end
    k = k - 11;
    for ll = 1:length(labels)
        [~,idx] = max(datamat_cells_norm{2}{c}, [], 2);
        [~,sortidx] = sort(idx);
        datamat_cellsSorted_norm{ll} = datamat_cells_norm{ll}{c}(sortidx, :);
    
        subplot(2,6,k)
        
        imagesc(datamat_cellsSorted_norm{ll});
        xticks([1,size(datamat_cellsSorted_norm{ll},2)])
        xticklabels([0,120])
        yticks([1,size(datamat_cellsSorted_norm{ll},1)])
        yticklabels([1,size(datamat_cellsSorted_norm{ll},1)])
        set(gca, 'FontSize', 14)
        xlabel('Position (cm)')
        if k == 1 || k == 7
            ylabel('Neuron')
        end
        k = k + 6;
    end
    k = k - 11;

end

cbar1 = colorbar('Position', [0.92 0.585 0.01 0.34]);
cbar1.Ticks = [0,1];
cbar1.FontSize = 12;
title(cbar1, 'Hz', 'FontSize', 12)

cbar2 = colorbar('Position', [0.92 0.11 0.01 0.34]);
cbar2.Ticks = [0,1];
cbar2.FontSize = 12;
title(cbar2, 'Hz', 'FontSize', 12)

annotation(figure(1),'textbox',...
    [0.0583824317086234 0.732585859488961 0.0372254944974622 0.0429292920261922],...
    'String',{'right'},'LineStyle','none','FontWeight','bold','FontSize',14);
annotation(figure(1),'textbox',...
    [0.0589180503481521 0.255313132216234 0.0297268337675385 0.0429292920261923],...
    'String',{'left'},'LineStyle','none','FontWeight','bold','FontSize',14);
annotation(figure(1),'textbox',...
    [0.215318693090519 0.943444445347551 0.0484734855923477 0.0429292920261922],...
    'String',{'2 ports'},'LineStyle','none', 'FontWeight','bold','FontSize',14, ...
    'Color',[0, 0.2471, 0.3608]);
annotation(figure(1),'textbox',...
    [0.487412961971076 0.938393940297045 0.0425816807331219 0.0429292920261922],...
    'String',{'1 port'},'LineStyle','none','FontWeight','bold','FontSize',14, ...
    'Color',[0.5843, 0.3176, 0.5882]);
annotation(figure(1),'textbox',...
    [0.750937332619175 0.942181819084922 0.0549009090751394 0.0429292920261922],...
    'String',{'walls off'},'LineStyle','none','FontWeight','bold','FontSize',14, ...
    'Color',[1, 0.6510, 0]);

saveas(gcf,[root filesep 'FiringMapAvg/PlaceCells',filesep,'heatmaps_norm.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/PlaceCells',filesep,'heatmaps_norm.fig'],'fig');


%% Correlate right vs left map 
titles = ["2 ports", "2 ports flipped", "1 port", "1 port flipped", "no walls", "no walls flipped"];

colors = [[0, 0.2471, 0.3608]; ...
    [0.3294, 0.4941, 0.6039]; ...
    [0.5843, 0.3176, 0.5882]; ...
    [0.7412, 0.5255, 0.7373]; ...
    [1, 0.6510, 0]; ...
    [1.0000, 0.7373, 0.3804]];

%% Cell by cell (row by row) 
for c = 1:3
    corrMapCell{c} = zeros(size(datamat_cells_norm{1}{c},1),1);
    corrMapCellFlip{c} = zeros(size(datamat_cells_norm{1}{c},1),1);
    for i = 1:size(datamat_cells_norm{1}{c},1) 
        corr = corrcoef(datamat_cells_norm{1}{c}(i,:)', datamat_cells_norm{2}{c}(i,:)','rows','complete');
        corrMapCell{c}(i) = corr(1,2);
        corrFlip = corrcoef(datamat_cells_norm{1}{c}(i,:)', flip(datamat_cells_norm{2}{c}(i,:))','rows','complete');
        corrMapCellFlip{c}(i) = corrFlip(1,2);
    end
end

% Statistical tests and plot
dataCombined = [corrMapCell{1}; corrMapCell{2}; corrMapCell{3}; ...
    corrMapCellFlip{1}; corrMapCellFlip{2}; corrMapCellFlip{3}];
group1 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{2}),1)*2; ...
    ones(length(corrMapCell{3}),1)*3; ones(length(corrMapCellFlip{1}),1); ...
    ones(length(corrMapCellFlip{2}),1)*2;ones(length(corrMapCellFlip{3}),1)*3];
group2 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{2}),1)*1; ...
    ones(length(corrMapCell{3}),1)*1; ones(length(corrMapCellFlip{1}),1)*2; ... 
    ones(length(corrMapCellFlip{2}),1)*2;ones(length(corrMapCellFlip{3}),1)*2];

directionStatsAvg.correlation.stats_all  = groupStats(dataCombined, [group1 group2], ...
    'color',colors,'plotType','boxplot','labelSummary',false,'sigStarTest','anova');
set(gcf,'Position',[300 100 700 600]);
ylabel('Corr. Coefficient', 'FontSize',14)
xticks([1,2,3.5,4.5,6,7])
xticklabels([titles(1), titles(2), titles(3), titles(4), titles(5), titles(6)])
yticks([-1,0,1])
set(gca, 'FontSize', 14)

saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_CellCell_new.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_CellCell_new.fig'],'fig');

% Separate by conditions (2 ports vs walls off, 2 ports vs 1 port)
titles1 = ["2 ports", "2 ports flipped", "no walls", "no walls flipped"];
titles2 = ["2 ports", "2 ports flipped", "1 port", "1 port flipped"];
colors1 = [[0, 0.2471, 0.3608]; ...
    [0.3294, 0.4941, 0.6039]; ...
    [1, 0.6510, 0]; ...
    [1.0000, 0.7373, 0.3804]];
colors2 = [[0, 0.2471, 0.3608]; ...
    [0.3294, 0.4941, 0.6039]; ...
    [0.5843, 0.3176, 0.5882]; ...
    [0.7412, 0.5255, 0.7373]];

% 2 ports vs walls off
dataCombined1 = [corrMapCell{1}; corrMapCell{3}; ...
    corrMapCellFlip{1}; corrMapCellFlip{3}];
group1 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{3}),1)*2; 
    ones(length(corrMapCellFlip{1}),1); ones(length(corrMapCellFlip{3}),1)*2];
group2 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{3}),1)*1; 
    ones(length(corrMapCellFlip{1}),1)*2; ones(length(corrMapCellFlip{3}),1)*2];

directionStatsAvg.correlation.stats_2portswallsoff = groupStats(dataCombined1, [group1 group2], ...
    'color',colors1,'plotType','boxplot','labelSummary',false,'sigStarTest','anova');
set(gcf,'Position',[2172,289,804,600]);
ylabel('Corr. Coefficient', 'FontSize',18)
xticks([1,2,3.5,4.5])
xticklabels([titles1(1), titles1(2), titles1(3), titles1(4)])
yticks([-1,0,1])
set(gca, 'FontSize', 18)
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_CellCell_2portsWallsOff.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_CellCell_2portsWallsOff.fig'],'fig');

% 2 ports vs 1 port
dataCombined2 = [corrMapCell{1}; corrMapCell{2}; ...
    corrMapCellFlip{1}; corrMapCellFlip{2}];
group1 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{2}),1)*2; 
    ones(length(corrMapCellFlip{1}),1); ones(length(corrMapCellFlip{2}),1)*2];
group2 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{2}),1)*1; 
    ones(length(corrMapCellFlip{1}),1)*2; ones(length(corrMapCellFlip{2}),1)*2];

directionStatsAvg.correlation.stats_2ports1port = groupStats(dataCombined2, [group1 group2], ...
    'color',colors2,'plotType','boxplot','labelSummary',false,'sigStarTest','anova');
set(gcf,'Position',[2103,289,873,600]);
ylabel('Corr. Coefficient', 'FontSize',18)
xticks([1,2,3.5,4.5])
xticklabels([titles2(1), titles2(2), titles2(3), titles2(4)])
yticks([-1,0,1])
set(gca, 'FontSize', 18)
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_CellCell_2ports1port.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_CellCell_2ports1port.fig'],'fig');

%% Position bin by bin (column by column) 
figure;
set(gcf,'Position',[300 100 800 200]);
k = 1;
j = 2;
for c = 1:3
    corrMapBin{c} = zeros(size(datamat_cells_norm{1}{c},2),1);
    corrMapBinFlip{c} = zeros(size(datamat_cells_norm{1}{c},2),1);
    for i = 1:size(datamat_cells_norm{1}{c},2) 
        corr = corrcoef(datamat_cells_norm{1}{c}(:,i), datamat_cells_norm{2}{c}(:,i),'rows','complete');
        corrFlip = corrcoef(datamat_cells_norm{1}{c}(:,i), flip(datamat_cells_norm{2}{c}(:,i)),'rows','complete');
        corrMapBin{c}(i) = corr(1,2);
        corrMapBinFlip{c}(i) = corrFlip(1,2);
    end
   
    subplot(3,2,k)
    imagesc(corrMapBin{c}');
    title(titles(c))
    subplot(3,2,j)
    imagesc(corrMapBinFlip{c}');
    title(titles(c+3))
    k = k + 2;
    j = j + 2;
end
caxis_min = min([corrMapBin{1}, corrMapBin{2}, corrMapBin{3}, corrMapBinFlip{1}, ...
    corrMapBinFlip{2}, corrMapBinFlip{3}], [], 'all');
caxis_max = max([corrMapBin{1}, corrMapBin{2}, corrMapBin{3}, corrMapBinFlip{1}, ...
    corrMapBinFlip{2}, corrMapBinFlip{3}], [], 'all');

k = 1;
j = 2;
for s = 1:3
    subplot(3,2,k);
    clim([caxis_min, caxis_max]);
    set(gca,'ytick',[])
    subplot(3,2,j);
    clim([caxis_min, caxis_max]);
    set(gca,'ytick',[])
    k = k + 2;
    j = j + 2;
end
subplot(3,2,6);
colorbar('Position', [0.93, 0.1, 0.02, 0.78], 'Location', 'eastoutside');

saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_BinBin_new.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_BinBin_new.fig'],'fig');


% Bin by bin line plots 
figure;
set(gcf,'Position',[300 100 1300 700]);
subplot(1,2,1)
hold on
for c = 1:3
    plot(1:length(corrMapBin{c}), corrMapBin{c}, 'Color', colors(c,:), 'DisplayName', ...
        titles(c), 'LineWidth', 2)
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility','off' )
legend('FontSize', 12)
xlim([0,length(corrMapBinFlip{1})])
ylim([caxis_min - 0.2, caxis_max + 0.2])
xlabel('Position bin')
ylabel('Corr. Coefficient')

subplot(1,2,2)
hold on
for c = 1:3
    plot(1:length(corrMapBinFlip{c}), corrMapBinFlip{c}, 'Color', colors(c,:), ...
        'DisplayName', titles(c+3), 'LineWidth', 2)
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility','off' )
legend('FontSize', 12);
xlim([0,length(corrMapBinFlip{1})])
ylim([caxis_min - 0.2, caxis_max + 0.2])
xlabel('Position bin')
ylabel('Corr. Coefficient')

saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_BinBin_lines_new.png'],'png');
saveas(gcf,[root filesep 'FiringMapAvg/correlations',filesep,'norm_BinBin_lines_new.fig'],'fig');

%% Save output
if saveMat
   save([root,filesep,'directionStatsAvg.cellinfo.mat'],'directionStatsAvg'); 
end
close all 
end
