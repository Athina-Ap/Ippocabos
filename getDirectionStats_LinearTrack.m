function [directionStats] = getDirectionStats_LinearTrack(varargin)

% (1) Extract cells with place fields and re-plot the place maps. 
% (2) Extract cells with fields in both directions. 
% (3) Perform correlation analyses when there are multiple session per day.

close all 

p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'thresFR',0.5,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
thresFR = p.Results.thresFR;

%% Load data
% Get session info
basename = bz_BasenameFromBasepath(basepath);
sessionInfo = load([basepath filesep basename '.sessionInfo.mat']);

% Get firing maps 
if ~isempty([basepath filesep basename 'firingMapsTrial.cellinfo.mat'])
    load([basepath filesep basename '.firingMapsTrial.cellinfo.mat']);
else
    disp('No rate maps in the this directory.')
    exit 
end

% Get place fields
if ~isempty([basepath filesep basename 'placeFields.cellinfo.mat'])
    load([basepath filesep basename '.placeFields.cellinfo.mat']);
else
    disp('No palce fields in this directory.')
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

% Get merge points
if ~isempty(dir([basepath filesep '*.MergePoints.events.mat'])) 
    disp('Loading MergePoints');
    file = dir([basepath filesep '*.MergePoints.events.mat']);
    load(file(1).name);
end

% Get spikes 
if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

if exist('MergePoints') & length(MergePoints.foldernames) > 1
    if ~isempty([basepath filesep basename 'firingMapsCond.cellinfo.mat'])
        load([basepath filesep basename '.firingMapsCond.cellinfo.mat']);
    end
end

%% Start analysis 
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    labels{1} = 'right';
    labels{2} = 'left';
    conditions{1} = '1port';
    conditions{2} = '2ports';

    for cond = 1:length(placeFieldStats)
        numcells = length(placeFieldStats{cond}.mapStats);
        directionStats{cond} = getDirectionStats(labels, numcells, placeFieldStats{cond});

        % Make plots
        labels_cond{1} = strcat([labels{1}, ' ', conditions{cond}]);
        labels_cond{2} = strcat([labels{2}, ' ', conditions{cond}]);
        [datamat_cells_norm{cond}] = plotPlaceCellMaps(directionStats{cond}, thresFR, firingMapsCond{cond}, ...
            tracking, behavTrials, spikes, labels, conditions{cond});
    end
else
    if strcmp(placeFieldStats.params.start, 'left')
        labels{1} = 'right';
        labels{2} = 'left';
    else
        labels{1} = 'left';
        labels{2} = 'right';
    end
    numcells = length(placeFieldStats.mapStats);
    directionStats = getDirectionStats(labels, numcells, placeFieldStats);

    % Make plots
    [datamat_cells_norm] = plotPlaceCellMaps(directionStats, thresFR, firingMaps, tracking, ...
        behavTrials, spikes, labels, 0);
end

%% Make extra plots for double sessions
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    placeCellsAll = unique([directionStats{1}.placeCells.right, directionStats{1}.placeCells.left, ...
        directionStats{2}.placeCells.right, directionStats{2}.placeCells.left]);
    figure;
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[1975,24,1445,968])
    j = 1;

    datamat_cells{1}{1} = [];
    datamat_cells{1}{2} = [];
    datamat_cells{2}{1} = [];
    datamat_cells{2}{2} = [];

    for cond = 1:length(conditions)
        for ll = 1:length(labels)
            for pf = 1:length(placeCellsAll)
                datamat_cells{cond}{ll} = [datamat_cells{cond}{ll}; firingMapsCond{cond}.rateMaps{placeCellsAll(pf)}{ll}];
            end
        end
    end
    
    % Normalize firing rate
    for cond = 1:length(conditions)
        datamat_cells_normAll{cond} = zeros(size(datamat_cells{cond}{1},1), 2*size(datamat_cells{cond}{1},2));
        for pf = 1:size(datamat_cells{cond}{1},1)
            datamat_cells_normAll{cond}(pf,:) = normalize([datamat_cells{cond}{1}(pf,:), datamat_cells{cond}{2}(pf,:)], 'range');
            norm_data = datamat_cells_normAll{cond}(pf,:);
            datamat_cells_norm{cond}{1}(pf,:) = norm_data(1:length(datamat_cells{cond}{1}(pf,:)));
            datamat_cells_norm{cond}{2}(pf,:) = norm_data(length(datamat_cells{cond}{1}(pf,:))+1:end);
        end
    end

    for cond = 1:length(conditions)
        for ll = 1:length(labels)
            [~,idx] = max(datamat_cells_norm{1}{1}, [], 2);
            [~,sortidx] = sort(idx);
            datamat_cellsSorted{cond}{ll} = datamat_cells_norm{cond}{ll}(sortidx, :);

            subplot(2,4,j)
%                 imagesc(zscore(datamat_cellsSorted.(labels{ll}){f},[],2));
            imagesc(datamat_cellsSorted{cond}{ll});
            xticks([1,size(datamat_cellsSorted{cond}{ll},2)])
            xticklabels([0,120])
            yticks([1,size(datamat_cellsSorted{cond}{ll},1)])
            yticklabels([1,size(datamat_cellsSorted{cond}{ll},1)])
            set(gca, 'FontSize', 14)
            xlabel('Position (cm)')
            if j == 1 || j == 5
                ylabel('Neuron')
            end
            j = j + 4;
        end
        j = j - 7;
        for ll = 1:length(labels)
            [~,idx] = max(datamat_cells_norm{1}{2}, [], 2);
            [~,sortidx] = sort(idx);
            datamat_cellsSorted{cond}{ll} = datamat_cells_norm{cond}{ll}(sortidx, :);
            subplot(2,4,j)
%                 imagesc(zscore(datamat_cellsSorted.(labels{ll}){f},[],2));
            imagesc(datamat_cellsSorted{cond}{ll});
            xticks([1,size(datamat_cellsSorted{cond}{ll},2)])
            xticklabels([0,120])
            yticks([1,size(datamat_cellsSorted{cond}{ll},1)])
            yticklabels([1,size(datamat_cellsSorted{cond}{ll},1)])
            set(gca, 'FontSize', 14)
            xlabel('Position (cm)')
            if j == 1 || j == 5
                ylabel('Neuron')
            end
            j = j + 4;
        end
        j = j - 7;
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
    [0.276058765129283 0.91749466059987 0.0648227848101265 0.0533707865168538],...
    'Color',[0.5843, 0.3176, 0.5882],...
    'String',{'1 port'},'LineStyle','none','FontWeight','bold', ...
    'FontSize',14,'FitBoxToText','off');

    annotation(figure(1),'textbox',...
    [0.688667352868175 0.912886526139844 0.0815089529482677 0.0533707865168541],...
    'Color',[0, 0.2471, 0.3608],...
    'String','2 ports','LineStyle','none','FontWeight','bold',...
    'FontSize',14,'FitBoxToText','off');

    annotation(figure(1),'textbox',...
    [0.0464528199962026 0.238299749280342 0.0455887024951468 0.0533707865168541],...
    'String','left','LineStyle','none','FontWeight','bold', ...
    'FontSize',14,'FitBoxToText','off');

    annotation(figure(1),'textbox',...
    [0.0443766954287285 0.719704707958028 0.0455887024951469 0.0533707865168541],...
    'String','right','LineStyle','none','FontWeight','bold',...
    'FontSize',14,'FitBoxToText','off');

    saveas(gcf,['FiringMapAvg/PlaceCells',filesep,'placeCells_2cond.png'],'png');
    saveas(gcf,['FiringMapAvg/PlaceCells',filesep,'placeCells_2cond.fig'],'fig');
    close(figure(1))
end


%% Correlation analysis 
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    
%     colors = [[0.1843, 0.2941, 0.4863]; ...
%         [0.4039, 0.3176, 0.5686]; ...
%         [0.9765, 0.3647, 0.4157]; ...
%         [1.0000, 0.4863, 0.2667]]; 
    colors = [[0.0000, 0.5451, 0.5451]; ...
        [0.3922, 0.5843, 0.9294]; ...
        [0.8039, 0.3608, 0.3608]; ...
        [1.0000, 0.2706, 0.0000]];

    titles = ["1 port", "2 ports", "1 port flipped", "2 ports flipped"];
    titles2 = {'same direction', 'same direction flipped', 'opposite direction', 'opposite direction flipped'};
    titles3 = ["same direction (R)", "same direction (L)", "opposite direction (RL)", ...
        "opposite direction (LR)", "same direction flipped (R)", "same direction flipped (L)", ...
        "opposite direction flipped (RL)", "opposite direction flipped (LR)"];
    
    %% Cell-cell correlation
    % Per direction within condition
    for cond = 1:length(conditions)
        corrMapCell{cond} = zeros(length(placeCellsAll),1);
        corrMapCellFlip{cond} = zeros(length(placeCellsAll),1);
    
        for i = 1:length(placeCellsAll)
            corr = corrcoef(datamat_cells_norm{cond}{1}(i,:)', datamat_cells_norm{cond}{2}(i,:)','rows','complete');
            corrMapCell{cond}(i) = corr(1,2);
            corrFlip = corrcoef(datamat_cells_norm{cond}{1}(i,:)', flip(datamat_cells_norm{cond}{2}(i,:))','rows','complete');
            corrMapCellFlip{cond}(i) = corrFlip(1,2);
        end
    end

    % Per condition within direction & per direction across conditions
    for ll = 1:length(labels)
        od = 3-ll;
        corrMapCellCond{ll} = zeros(length(placeCellsAll),1);
        corrMapCellFlipCond{ll} = zeros(length(placeCellsAll),1);
        corrMapCellDir{ll} = zeros(length(placeCellsAll),1);
        corrMapCellFlipDir{ll} = zeros(length(placeCellsAll),1);

        for i = 1:length(placeCellsAll)
            corr = corrcoef(datamat_cells_norm{1}{ll}(i,:)', datamat_cells_norm{2}{ll}(i,:)','rows','complete');
            corrMapCellCond{ll}(i) = corr(1,2);
            corrFlip = corrcoef(datamat_cells_norm{1}{ll}(i,:)', flip(datamat_cells_norm{2}{ll}(i,:))','rows','complete');
            corrMapCellFlipCond{ll}(i) = corrFlip(1,2);

            corr = corrcoef(datamat_cells_norm{1}{ll}(i,:)', datamat_cells_norm{2}{od}(i,:)','rows','complete');
            corrMapCellDir{ll}(i) = corr(1,2);
            corrFlip = corrcoef(datamat_cells_norm{1}{ll}(i,:)', flip(datamat_cells_norm{2}{od}(i,:))','rows','complete');
            corrMapCellFlipDir{ll}(i) = corrFlip(1,2);
        end
    end
    
    % Statistical tests and plot
    % Per direction within condition
    dataCombined = [corrMapCell{1}; corrMapCell{2}; ...
        corrMapCellFlip{1}; corrMapCellFlip{2}; ];
    group1 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{2}),1)*2; ...
        ones(length(corrMapCellFlip{1}),1); ones(length(corrMapCellFlip{2}),1)*2;]; % group by condition
    group2 = [ones(length(corrMapCell{1}),1); ones(length(corrMapCell{2}),1)*1; ...
        ones(length(corrMapCellFlip{1}),1)*2; ones(length(corrMapCellFlip{2}),1)*2]; % group by flipped or normal
    
    directionStats{3}.CellCellCorr.stats = groupStats(dataCombined, [group1 group2], ...
        'color',colors,'plotType','boxplot','labelSummary',false,'sigStarTest','anova');
    set(gcf,'Position',[300 100 700 600]);
    ylabel('Corr. Coefficient')
    xticks([1,2,3.5,4.5,6,7])
    xticklabels([titles(1), titles(2), titles(3), titles(4)])
    
    if ~isdir('FiringMapAvg\correlations')
        mkdir('FiringMapAvg\correlations')
    end
    saveas(gcf,[basepath filesep 'FiringMapAvg\correlations',filesep,'CellCell.png'],'png');
    saveas(gcf,[basepath filesep 'FiringMapAvg\correlations',filesep,'CellCell.fig'],'fig');

    % Per condition within direction & per direction across conditions
    data1 = [corrMapCellCond{1}; corrMapCellCond{2}];
    data2 = [corrMapCellDir{1}; corrMapCellDir{2}];
    data3 = [corrMapCellFlipCond{1}; corrMapCellFlipCond{2}];
    data4 = [corrMapCellFlipDir{1}; corrMapCellFlipDir{2}];
    dataCombined = [data1; data2; data3; data4];

    group1 = [ones(length(data1),1); ones(length(data2),1)*2; ...
        ones(length(data3),1); ones(length(data4),1)*2]; % group by direction 

    group2 = [ones(length(data1),1); ones(length(data2),1); ...
        ones(length(data3),1)*2; ones(length(data4),1)*2]; % group by flipped or normal
        
    directionStats{4}.CellCellCorr.stats = groupStats(dataCombined, [group1 group2], ...
        'color',colors,'plotType','boxplot','labelSummary',false,'sigStarTest','anova');
    set(gcf,'Position',[0.3922, 0.5843, 0.9294]);
    ylabel('Corr. Coefficient')
    xticks([1,2,3.5,4.5,6,7])
    yticks([-1,0,1])
    set(gca, 'FontSize', 18)
    xticklabels([titles2(1), titles2(2), titles2(3), titles2(4)])
    
    saveas(gcf,[basepath filesep 'FiringMapAvg\correlations',filesep,'CellCell_conditions.png'],'png');
    saveas(gcf,[basepath filesep 'FiringMapAvg\correlations',filesep,'CellCell_conditions.fig'],'fig');

    %% Bin-bin correlations
    numbins = size(datamat_cells_norm{cond}{1},2);

    % Per direction within condition 
    for cond = 1:length(conditions)
        corrMapBin{cond} = zeros(numbins,1);
        corrMapBinFlip{cond} = zeros(numbins,1);
        for i = 1:numbins
            corr = corrcoef(datamat_cells_norm{cond}{1}(:,i), datamat_cells_norm{cond}{2}(:,i),'rows','complete');
            corrFlip = corrcoef(datamat_cells_norm{cond}{1}(:,i), flip(datamat_cells_norm{cond}{2}(:,i)),'rows','complete');
            corrMapBin{cond}(i) = corr(1,2);
            corrMapBinFlip{cond}(i) = corrFlip(1,2);
        end
    end
    
    % Per condition within direction & per direction across conditions
    for ll = 1:length(labels)
        od = 3-ll;
        corrMapBinCond{ll} = zeros(numbins,1);
        corrMapBinFlipCond{ll} = zeros(numbins,1);
        corrMapBinDir{ll} = zeros(numbins,1);
        corrMapBinFlipDir{ll} = zeros(numbins,1);

        for i = 1:numbins
            corr = corrcoef(datamat_cells_norm{1}{ll}(:,i)', datamat_cells_norm{2}{ll}(:,i)','rows','complete');
            corrMapBinCond{ll}(i) = corr(1,2);
            corrFlip = corrcoef(datamat_cells_norm{1}{ll}(:,i)', flip(datamat_cells_norm{2}{ll}(:,i))','rows','complete');
            corrMapBinFlipCond{ll}(i) = corrFlip(1,2);

            corr = corrcoef(datamat_cells_norm{1}{ll}(:,i)', datamat_cells_norm{2}{od}(:,i)','rows','complete');
            corrMapBinDir{ll}(i) = corr(1,2);
            corrFlip = corrcoef(datamat_cells_norm{1}{ll}(:,i)', flip(datamat_cells_norm{2}{od}(:,i))','rows','complete');
            corrMapBinFlipDir{ll}(i) = corrFlip(1,2);
        end
    end

    % Heatmaps: per direction within condition 
    figure;
    set(gcf,'Position',[300 100 800 200]);
    k = 1;
    j = 2;
    for cond = 1:length(conditions)
        subplot(2,2,k)
        imagesc(corrMapBin{cond}');
        title(titles(cond))
        subplot(2,2,j)
        imagesc(corrMapBinFlip{cond}');
        title(titles(cond+2))
        k = k + 2;
        j = j + 2;
    end
    caxis_min = min([corrMapBin{1}, corrMapBin{2}, corrMapBinFlip{1}, corrMapBinFlip{2}], [], 'all');
    caxis_max = max([corrMapBin{1}, corrMapBin{2}, corrMapBinFlip{1}, corrMapBinFlip{2}], [], 'all');

    k = 1;
    j = 2;
    for s = 1:2
        subplot(2,2,k);
        clim([caxis_min, caxis_max]);
        set(gca,'ytick',[])
        subplot(2,2,j);
        clim([caxis_min, caxis_max]);
        set(gca,'ytick',[])
        k = k + 2;
        j = j + 2;
    end
    subplot(2,2,4);
    colorbar('Position', [0.93, 0.1, 0.02, 0.78], 'Location', 'eastoutside');

    saveas(gcf,[basepath filesep 'FiringMapAvg/correlations',filesep,'BinBin_heatmap.png'],'png');
    saveas(gcf,[basepath filesep 'FiringMapAvg/correlations',filesep,'BinBin_heatmap.fig'],'fig');

    % Line plots
    % Per direction within condition 
    figure;
    set(gcf,'Position',[300 100 1300 700]);
    subplot(1,2,1)
    hold on
    for cond = 1:length(conditions)
        plot(1:length(corrMapBin{cond}), corrMapBin{cond}, 'Color', colors(cond,:), 'DisplayName', ...
            titles(cond), 'LineWidth', 2)
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility','off' )
    legend('FontSize', 12)
    xlim([0,length(corrMapBinFlip{1})])
    ylim([caxis_min - 0.2, caxis_max + 0.2])
    xlabel('Position bin')
    ylabel('Corr. Coefficient')
    
    subplot(1,2,2)
    hold on
    for cond = 1:length(conditions)
        plot(1:length(corrMapBinFlip{cond}), corrMapBinFlip{cond}, 'Color', colors(cond,:), ...
            'DisplayName', titles(cond+2), 'LineWidth', 2)
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility','off' )
    legend('FontSize', 12);
    xlim([0,length(corrMapBinFlip{1})])
    ylim([caxis_min - 0.2, caxis_max + 0.2])
    xlabel('Position bin')
    ylabel('Corr. Coefficient')
    
    saveas(gcf,[basepath filesep 'FiringMapAvg/correlations',filesep,'BinBin_lines.png'],'png');
    saveas(gcf,[basepath filesep 'FiringMapAvg/correlations',filesep,'BinBin_lines.fig'],'fig');


    % Per condition within direction & per direction across conditions
    figure;
    set(gcf,'Position',[300 100 1300 700]);
    subplot(1,2,1)
    hold on
    j = 1;
    for ll = 1:length(labels)
        plot(1:numbins, corrMapBinCond{ll}, 'Color', colors(j,:), 'DisplayName', ...
            titles3(j), 'LineWidth', 2)
        plot(1:numbins, corrMapBinDir{ll}, 'Color', colors(j+2,:), 'DisplayName', ...
            titles3(j+2), 'LineWidth', 2)
        j = j + 1;
    end

    yline(0, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility','off' )
    legend('FontSize', 12)
    xlim([0,numbins])
    ylim([caxis_min - 0.2, caxis_max + 0.2])
    xlabel('Position bin')
    ylabel('Corr. Coefficient')
    
    subplot(1,2,2)
    hold on
    j = 1;
    for ll = 1:length(labels)
        plot(1:numbins, corrMapBinFlipCond{ll}, 'Color', colors(j,:), ...
            'DisplayName', titles3(j+4), 'LineWidth', 2)
        plot(1:numbins, corrMapBinFlipDir{ll}, 'Color', colors(j+2,:), ...
            'DisplayName', titles3(j+6), 'LineWidth', 2)
        j = j + 1;
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility','off' )
    legend('FontSize', 12);
    xlim([0,numbins])
    ylim([caxis_min - 0.2, caxis_max + 0.2])
    xlabel('Position bin')
    ylabel('Corr. Coefficient')

    saveas(gcf,[basepath filesep 'FiringMapAvg/correlations',filesep,'BinBin_lines_conditions.png'],'png');
    saveas(gcf,[basepath filesep 'FiringMapAvg/correlations',filesep,'BinBin_lines_conditions.fig'],'fig');
    close all
end

%% Save the output 
if saveMat
   save([basepath,filesep,basename '.directionStats.cellinfo.mat'],'directionStats'); 
end
end
