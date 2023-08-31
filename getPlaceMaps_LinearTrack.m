function getPlaceMaps_LinearTrack(varargin)

close all

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'frThres',true,@islogical);
addParameter(p,'thres',0.5,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
plotfig = p.Results.plotfig;
frThres = p.Results.frThres;
thres = p.Results.thres;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.MergePoints.events.mat'])) 
    disp('Loading MergePoints');
    file = dir([basepath filesep '*.MergePoints.events.mat']);
    load(file(1).name);
end
    
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end
  
if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

% Correct too frequent solenoid breaks
behavTrials.timestampsNew = behavTrials.timestamps(:)';
delays = diff(behavTrials.timestampsNew); 
errorsIdx = find(delays < 2);
behavTrials.timestampsNew(errorsIdx) = [];

%% Compute place fields
% Compute two sets of maps. Left running, right running. 
fprintf('Computing place fields.\n');

%% Assign spike position to each spike
if ~isempty(dir([basepath filesep '*.spikeData.cellinfo.mat']))
    file = dir([basepath filesep '*.spikeData.cellinfo.mat']);
    load(file.name);
else
    for unit = 1:length(spikes.UID)
        [idx] = InIntervals(spikes.times{unit},[tracking.timestamps(1) tracking.timestamps(end)]);
        tsBehav = spikes.times{unit}(idx);
        for tt = 1:length(tsBehav)
            [~,closestIndex] = min(abs(tracking.timestamps-tsBehav(tt)));
            spikeData.posIdx{unit}(tt) = closestIndex;
        end
        spikeData.pos{unit} = tracking.position.y(spikeData.posIdx{unit});
    end
    save([sessionInfo.FileName '.spikeData.cellinfo.mat'],'spikeData'); 
end

%% Bin positions in left vs right running directions
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    if length(MergePoints.foldernames) == 2
        positions.right = [];
        positions.left = [];
        idxR_trackedAll = [];
        idxL_trackedAll = [];

        mergeidx = find(behavTrials.timestamps(:,1) > MergePoints.timestamps(1,2), 1,'first');

        for f = 1:length(MergePoints.foldernames)
            if f == 1
                sessLen = 1:mergeidx;
                [positions_single, idxL_tracked{f}, idxR_tracked{f}] = getTrackedPositions(behavTrials, ...
                    tracking, f, sessLen);
                positions.right = [positions.right, positions_single.right];
                positions.left = [positions.left, positions_single.left];
                idxR_trackedAll = [idxR_trackedAll, idxR_tracked{f}];
                idxL_trackedAll = [idxL_trackedAll, idxL_tracked{f}];
                
                positions.allTrialsCond{f}{1} = [tracking.timestamps(idxR_tracked{f}) tracking.position.y(idxR_tracked{f})];
                positions.allTrialsCond{f}{2} = [tracking.timestamps(idxL_tracked{f}) tracking.position.y(idxL_tracked{f})];
            
            elseif f == 2
                sessLen = mergeidx:size(behavTrials.timestamps,1);
                [positions_single, idxL_tracked{f}, idxR_tracked{f}] = getTrackedPositions(behavTrials, ...
                    tracking, f, sessLen);
                positions.right = [positions.right, positions_single.right];
                positions.left = [positions.left, positions_single.left];
                idxR_trackedAll = [idxR_trackedAll, idxR_tracked{f}];
                idxL_trackedAll = [idxL_trackedAll, idxL_tracked{f}];
                
                positions.allTrialsCond{f}{1} = [tracking.timestamps(idxR_tracked{f}) tracking.position.y(idxR_tracked{f})];
                positions.allTrialsCond{f}{2} = [tracking.timestamps(idxL_tracked{f}) tracking.position.y(idxL_tracked{f})];
            end
        end
        positions.allTrials{1} = [tracking.timestamps(idxR_trackedAll) tracking.position.y(idxR_trackedAll)];
        positions.allTrials{2} = [tracking.timestamps(idxL_trackedAll) tracking.position.y(idxL_trackedAll)];
    end                
else
    sessLen = 1:size(behavTrials.timestamps,1);
    [positions, idxL_tracked, idxR_tracked] = getTrackedPositions(behavTrials, ...
                    tracking, 0, sessLen);
    if strcmp(behavTrials.start, 'left')
        positions.allTrials{1} = [tracking.timestamps(idxR_tracked) tracking.position.y(idxR_tracked)]; 
        positions.allTrials{2} = [tracking.timestamps(idxL_tracked) tracking.position.y(idxL_tracked)];
    else
        positions.allTrials{1} = [tracking.timestamps(idxL_tracked) tracking.position.y(idxL_tracked)];
        positions.allTrials{2} = [tracking.timestamps(idxR_tracked) tracking.position.y(idxR_tracked)];
    end
end


%% Save all trials together to determine if the cell has a place field
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    firingMapsCond{1} = bz_firingMapAvg_IZ(positions.allTrialsCond{1},spikes, ...
        'minTime',0.0001,'plotFig',false,'saveMat',false);
    firingMapsCond{2} = bz_firingMapAvg_IZ(positions.allTrialsCond{2},spikes, ...
        'minTime',0.0001,'plotFig',false,'saveMat',false);
    firingMapsCond{1}.start = behavTrials.start(1);
    firingMapsCond{2}.start = behavTrials.start(2); 
end
% firingMaps = bz_firingMapAvg_IZ(positions.allTrials,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);
% firingMaps.start = behavTrials.start;

firingMaps.right = bz_firingMapAvg_IZ(positions.right,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);
firingMaps.left = bz_firingMapAvg_IZ(positions.left,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);

save([sessionInfo.FileName '.firingMapsTrial.cellinfo.mat'],'firingMaps');
save([sessionInfo.FileName '.firingMapsCond.cellinfo.mat'], 'firingMapsCond');

%% Plot place maps
if exist('MergePoints') & length(MergePoints.foldernames) > 1
    labels = {'right', 'left'}; 
else
    if strcmp(behavTrials.start, 'left')
        labels = {'right', 'left'}; 
    else
        labels = {'left', 'right'};
    end
end

numtrials = length(firingMaps.right.rateMaps{1,1});
numcells = length(firingMaps.right.rateMaps);

clear idx
if plotfig
    
    %% (1) Plot all trials for each cell 
    if ~isfolder('FiringMap/CellTrial')
        mkdir('FiringMap/CellTrial')
    end
   
    for pf = 1:numcells
         figure
         set(gcf,'Renderer','painters')
         set(gcf,'Position',[2200 200 1185 712])
         datamat{1} = [];
         datamat{2} = [];
         for ll = 1:length(labels)
            for kk = 1:numtrials
                datamat{ll} = [datamat{ll}; firingMaps.(labels{ll}).rateMaps{pf}{kk}];
            end
            subplot(1,2,ll) 

            if ~isempty(datamat{ll})
                h = imagesc(datamat{ll});
                set(h, 'AlphaData', ~isnan(datamat{ll}))
                title(labels{ll})
                colorbar;
                
            else
                axis off
            end
         end                  
    
        saveas(gcf,['FiringMap/CellTrial',filesep,'cell_' num2str(pf) '.png'],'png');
        saveas(gcf,['FiringMap/CellTrial',filesep,'cell_' num2str(pf) '.fig'],'fig');
        close all;
    end

    %% (2) Plot the firing of each cell along the trajectory 
    getTrajectoryFiring_LinearTrack(1:numcells);

    %% (3) Plot the average across trials for each cell (per session)
    if ~isfolder('FiringMap/CellTrialAvg')
        mkdir('FiringMap/CellTrialAvg')
    end
    
    if exist('MergePoints') & length(MergePoints.foldernames) > 1
        datamat_mean.right{1} = cell(length(1:mergeidx),1);
        datamat_mean.left{1} = cell(length(1:mergeidx),1);
        datamat_mean.right{2} = cell(length(mergeidx:size(behavTrials.timestamps,1)),1);
        datamat_mean.left{2} = cell(length(mergeidx:size(behavTrials.timestamps,1)),1);

        for pf = 1:numcells  
            figure
            set(gcf,'Renderer','painters')
            set(gcf,'Position',[2200 200 1185 712])

            datamat_split.right{1} = [];
            datamat_split.right{2} = [];
            datamat_split.left{1} = [];
            datamat_split.left{2} = [];
            
            j = 1;
            for f = 1:length(MergePoints.foldernames)
                if f == 1
                    for ll = 1:length(labels)
                       for kk = 1:mergeidx
                           datamat_split.(labels{ll}){f} = [datamat_split.(labels{ll}){f}; firingMaps.(labels{ll}).rateMaps{pf}{kk}];
                       end
                       datamat_mean.(labels{ll}){f}{pf} = mean(datamat_split.(labels{ll}){f}, 1); 

                       subplot(2,2,j)
                       if ~isempty(datamat_mean.(labels{ll}){f})
                           h = imagesc(datamat_mean.(labels{ll}){f}{pf});
                           set(h, 'AlphaData', ~isnan(datamat_mean.(labels{ll}){f}{pf}))
                           set(gca, 'YTick', []);
                           set(gca, 'YTickLabel', []);
                           title(labels{ll})
                           colorbar;   
                       else
                           axis off
                       end
                       j = j + 1;
                    end
                elseif f == 2
                    for ll = 1:length(labels)
                       for kk = mergeidx:size(behavTrials.timestamps,1)-1
                           datamat_split.(labels{ll}){f} = [datamat_split.(labels{ll}){f}; firingMaps.(labels{ll}).rateMaps{pf}{kk}];
                       end
                       datamat_mean.(labels{ll}){f}{pf} = mean(datamat_split.(labels{ll}){f}, 1); 

                       subplot(2,2,j)
                       if ~isempty(datamat_mean.(labels{ll}){f})
                           h = imagesc(datamat_mean.(labels{ll}){f}{pf});
                           set(h, 'AlphaData', ~isnan(datamat_mean.(labels{ll}){f}{pf}))
                           set(gca, 'YTick', []);
                           set(gca, 'YTickLabel', []);
                           title(labels{ll})
                           colorbar;   
                       else
                           axis off
                       end
                       j = j + 1;
                    end
                end
            end       
           saveas(gcf,['FiringMap/CellTrialAvg',filesep,'cell_' num2str(pf) '.png'],'png');
           saveas(gcf,['FiringMap/CellTrialAvg',filesep,'cell_' num2str(pf) '.fig'],'fig');
           close all;
        end
    else
        datamat_mean{1} = cell(length(firingMaps.right.rateMaps),1);
        datamat_mean{2} = cell(length(firingMaps.right.rateMaps),1);

        for pf = 1:numcells
            figure
            set(gcf,'Renderer','painters')
            set(gcf,'Position',[2200 200 1185 200])
             
            datamat{1} = [];
            datamat{2} = [];
    
            for ll = 1:length(labels)
               for kk = 1:numtrials
                   datamat{ll} = [datamat{ll}; firingMaps.(labels{ll}).rateMaps{pf}{kk}];
               end
               datamat_mean{ll}{pf} = mean(datamat{ll}, 1); 
    
               subplot(1,2,ll) 
    
               if ~isempty(datamat_mean{ll})
                   h = imagesc(datamat_mean{ll}{pf});
                   set(h, 'AlphaData', ~isnan(datamat_mean{ll}{pf}))
                   set(gca, 'YTick', []);
                   set(gca, 'YTickLabel', []);
                   title(labels{ll})
                   colorbar;       
               else
                   axis off
               end
            end                  
           saveas(gcf,['FiringMap/CellTrialAvg',filesep,'cell_' num2str(pf) '.png'],'png');
           saveas(gcf,['FiringMap/CellTrialAvg',filesep,'cell_' num2str(pf) '.fig'],'fig');
           close all;
        end
    end

    %% (4) Plot the average across trials for all cells
    % Sort cells based on where max firing occurs (trials mean)
    if ~isfolder('FiringMapAvg')
        mkdir('FiringMapAvg');
    end

    if exist('MergePoints') & length(MergePoints.foldernames) > 1
        if length(MergePoints.foldernames) == 2
            datamat_cells.right{1} = [];
            datamat_cells.right{2} = [];
            datamat_cells.left{1} = [];
            datamat_cells.left{2} = [];

            for f = 1:length(MergePoints.foldernames)
                if f == 1
                    for ll = 1:length(labels)
                        for pf = 1:numcells
                            datamat_cells.(labels{ll}){f} = [datamat_cells.(labels{ll}){f}; datamat_mean.(labels{ll}){f}{pf}];
                        end
                    end
                elseif f == 2
                    for ll = 1:length(labels)
                        for pf = 1:numcells
                            datamat_cells.(labels{ll}){f} = [datamat_cells.(labels{ll}){f}; datamat_mean.(labels{ll}){f}{pf}];
                        end
                    end
                end
            end
        end
    else
        datamat_cells{1} = [];
        datamat_cells{2} = [];
        for ll = 1:length(labels)
            for pf = 1:numcells
                datamat_cells{ll} = [datamat_cells{ll}; datamat_mean{ll}{pf}];
            end
        end
    end

    figure;
    if exist('MergePoints') & length(MergePoints.foldernames) > 1
        j = 1;
        for f = 1:length(MergePoints.foldernames)
            for ll = 1:length(labels)
                [~,idx] = max(datamat_cells.right{1}, [], 2);
                [~,sortidx] = sort(idx);
                datamat_cellsSorted.(labels{ll}){f} = datamat_cells.(labels{ll}){f}(sortidx, :);

                subplot(4,2,j)
                set(gcf,'Renderer','painters')
                set(gcf,'Position',[2200 200 1185 712])
                imagesc(zscore(datamat_cellsSorted.(labels{ll}){f},[],2));
                colorbar
                xlabel('Position bin')
                ylabel('Neuron')
                title(labels{ll})
                j = j + 1;
            end
            for ll = 1:length(labels)
                [~,idx] = max(datamat_cells.left{1}, [], 2);
                [~,sortidx] = sort(idx);
                datamat_cellsSorted.(labels{ll}){f} = datamat_cells.(labels{ll}){f}(sortidx, :);

                subplot(4,2,j)
                set(gcf,'Renderer','painters')
                set(gcf,'Position',[2200 200 1185 712])
                imagesc(zscore(datamat_cellsSorted.(labels{ll}){f},[],2));
                colorbar
                xlabel('Position bin')
                ylabel('Neuron')
                title(labels{ll})
                j = j + 1;
            end
        end
        annotation(figure(1),'textbox',...
        [0.0296919831223629 0.703651685393258 0.0648227848101266 0.0533707865168539],...
        'Color',[1 0.411764705882353 0.16078431372549],...
        'String',{'1 port'},'LineStyle','none','FontWeight','bold', ...
        'FontSize',14,'FitBoxToText','off');

        annotation(figure(1),'textbox',...
        [0.030535864978903 0.26825842696629 0.0648227848101265 0.053370786516854],...
        'Color',[0 0.447058823529412 0.741176470588235],...
        'String','2 ports','LineStyle','none','FontWeight','bold',...
        'FontSize',14,'FitBoxToText','off');
    else
        for ll = 1:length(labels)
            [~,idx] = max(datamat_cells{1}, [], 2);
            [~,sortidx] = sort(idx);
            datamat_cellsSorted{ll} = datamat_cells{ll}(sortidx, :);
    
            subplot(2,2,ll+2)
            set(gcf,'Renderer','painters')
            set(gcf,'Position',[2200 200 1185 712])
            imagesc(zscore(datamat_cellsSorted{ll},[],2));
            colorbar
            xlabel('Position bin')
            ylabel('Neuron')
            title(labels{ll})
        end
        for ll = 1:length(labels)
            [~,idx] = max(datamat_cells{2}, [], 2);
            [~,sortidx] = sort(idx);
            datamat_cellsSorted{ll} = datamat_cells{ll}(sortidx, :);
    
            subplot(2,2,ll+2)
            set(gcf,'Renderer','painters')
            set(gcf,'Position',[2200 200 1185 712])
            imagesc(zscore(datamat_cellsSorted{ll},[],2));
            colorbar
            xlabel('Position bin')
            ylabel('Neuron')
            title(labels{ll})
        end
    end
    
    saveas(gcf,['FiringMapAvg',filesep,'all_cells.png'],'png');
    saveas(gcf,['FiringMapAvg',filesep,'all_cells.fig'],'fig');
    close all;
end     

%% Apply a firing rate threshold
% 
% if frThres
%     dt = mean(diff(tracking.timestamps));
%     win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
%     spkData = bz_SpktToSpkmat(spikes,'dt',dt,'win',win);
%     
%     spkMat = spkData.data';
%     timestamps = spkData.timestamps';
% end
% 
% fr = 1:length(spkMat(:,1));
% for i = 1:length(spkMat(:,1))
%     fr(i) = mean(spkMat(i,:)) * 1/dt;
% end
% LowFiringCells = find(fr < thres);
% 
% % Plot firing maps
% if plotfig
%     datamat_cells{1} = [];
%     datamat_cells{2} = [];
%     
%     for ll = 1:length(labels)
%         for pf = 1:numcells
%             if ismember(pf, LowFiringCells)
%                 continue
%             end
%             datamat_cells{ll} = [datamat_cells{ll}; datamat_mean{ll}{pf}];
%         end
%     end
%     
%     figure;
%     for ll = 1:length(labels)
%         [~,idx] = max(datamat_cells{1}, [], 2);
%         [~,sortidx] = sort(idx);
%         datamat_cellsSorted{ll} = datamat_cells{ll}(sortidx, :);
%     
%         subplot(2,2,ll)
%         set(gcf,'Renderer','painters')
%         set(gcf,'Position',[2200 200 1185 712])
%         imagesc(zscore(datamat_cellsSorted{ll},[],2));
%         colorbar
%         xlabel('Position bin')
%         ylabel('Neuron')
%         title(labels{ll})
%     end
%     for ll = 1:length(labels)
%         [~,idx] = max(datamat_cells{2}, [], 2);
%         [~,sortidx] = sort(idx);
%         datamat_cellsSorted{ll} = datamat_cells{ll}(sortidx, :);
%     
%         subplot(2,2,ll+2)
%         set(gcf,'Renderer','painters')
%         set(gcf,'Position',[2200 200 1185 712])
%         imagesc(zscore(datamat_cellsSorted{ll},[],2));
%         colorbar
%         xlabel('Position bin')
%         ylabel('Neuron')
%         title(labels{ll})
%     end
%     
%     saveas(gcf,['FiringMapAvg',filesep,'all_cells_filt.png'],'png');
%     saveas(gcf,['FiringMapAvg',filesep,'all_cells_filt.fig'],'fig');
%     close all;
% end

end

