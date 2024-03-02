% This function calculates the percentage of cells for the chosen session
% that are tuned to each of the variables that we have chosen to test with
% the PGAM. 
% The function also identifies the cells that are tuned to licks and plots
% their PSTH. 

function [resultsPGAM] = postprocessPGAM_newVars2(varargin)

close all 
format long g

p = inputParser;
addParameter(p,'results','I:\April 2023 GAM fits\results\new_vars2',@isstr);
addParameter(p,'data','I:\Videos',@isstr);
addParameter(p,'postprocess','I:\April 2023 GAM fits\postprocess',@isstr);
% addParameter(p,'mapDir','Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ44_220830_sess7\Maps',@isstr);
parse(p,varargin{:});
results = p.Results.results;
data = p.Results.data;
postprocess = p.Results.postprocess;
% mapDir = p.Results.mapDir;

if ~isfolder(postprocess)
    mkdir(postprocess)
end

% sessions = {'IZ39\Final\IZ39_220622_sess8', 'IZ39\Final\IZ39_220624_sess10', 'IZ39\Final\IZ39_220629_sess12',...
%     'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
%     'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...  
%     'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
%     'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
%     'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
%     'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
% 'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
% 'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
%     'IZ44\Final\IZ44_220827_sess4','IZ44\Final\IZ44_220828_sess5',...
%     'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
%  'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
%  'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15'};
sessions = {'IZ47\Final\IZ47_230707_sess24', 'IZ47\Final\IZ47_230710_sess25', ...
    'IZ47\Final\230712_sess27', 'IZ48\Final\IZ48_230705_sess22'};

for s = 1:length(sessions)
    session = extractAfter(sessions{s}, 'Final\');
    dataDir = fullfile(data, sessions{s});
    resultsDir = fullfile(results, session);
    postprocessDir = fullfile(postprocess, session);
    
    % Load spiking data
    if ~isempty(dir([dataDir filesep '*.spikes.cellinfo.mat']))
        disp('Spikes already detected! Loading file.');
        file = dir([dataDir filesep '*.spikes.cellinfo.mat']);
        load(fullfile(dataDir, file.name));
    end
    
    % Load behavior data
    if ~isempty(dir([dataDir filesep '*TrialBehavior.Behavior.mat'])) 
        disp('Behavior already detected! Loading file.');
        file = dir([dataDir filesep '*TrialBehavior.Behavior.mat']);
        load(fullfile(dataDir, file.name));
    end
    
    % Load event data 
    if ~isempty(dir([dataDir filesep '*.sessionDataPGAM.mat']))
        disp('Variables already detected! Loading file.');
        file = dir([dataDir filesep '*.sessionDataPGAM.mat']);
        load(fullfile(dataDir, file.name));
    end
    
    % Load cell metrics
    if ~isempty(dir([dataDir filesep '*.cell_metrics.cellinfo.mat']))
        disp('Cell metrics already detected! Loading file.');
        file = dir([dataDir filesep '*.cell_metrics.cellinfo.mat']);
        load(fullfile(dataDir, file.name));
    end
    
    % Load digitalIn
    if ~isempty(dir([dataDir filesep '*.DigitalIn.events.mat']))
        disp('Digital input already detected! Loading file.');
        file = dir([dataDir filesep '*.DigitalIn.events.mat']);
        load(fullfile(dataDir, file.name));
    end
    
    % FIXME: numNeurons is different to number of files
    % Only select pyramidal neurons 
    resultsPGAM = {};
    
    numNeurons = spikes.numcells;
    pyramidal = find(contains(cell_metrics.putativeCellType, 'Pyramidal Cell'));
    
    resultsPGAM.numNeurons = numNeurons;
    resultsPGAM.pyramidal = pyramidal';
    
    % Load PGAM results
    filelist = dir(fullfile(resultsDir, '*k-fold.csv'));
    filename = {filelist.name}';
    
    %% Find the cells tuned to licks
    % Note here the indices refer to the pyramidal cells and not the cell ID.
    tuned_licks = zeros(length(pyramidal),1);
    tuned_relDistStop = zeros(length(pyramidal),1);
    
    for i = 1:length(pyramidal)
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(pyramidal(i)), '_fit_k-fold.csv')))
            
            results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
                string(pyramidal(i)), '_fit_k-fold.csv')));
        else
            continue
        end
    
        % Variables licks is in the last row (11) and the pval is in column 10
        if table2array(results(9, 10)) < 0.001
            tuned_licks(i, 1) = 1; 
        end
        if table2array(results(3, 10)) < 0.001
            tuned_relDistStop(i, 1) = 1; 
        end
    end
    resultsPGAM.tunedLicks = tuned_licks;
    
    %% Plot heatmap of cells tuned to licks - comment option 1/2 and change win
    dt = 0.033;
    win = 60;  % number of bins that capture 2 seconds
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine values needed to calculate the PSTH
    win2 = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
    spkData = bz_SpktToSpkmat(spikes,'dt',dt,'win',win2);
    
    spkMat = spkData.data';
    spk_timestamps = spkData.timestamps';
    
    lickport = [1 2 3 4 5 6 7];
    for ll = 1:7
        event.times{ll} = digitalIn.timestampsOn{lickport(ll)+2}; 
        event.UID(ll) = ll;
    end
    
    eventMat = bz_SpktToSpkmat(event,'dt',dt,'win',win2);
    event_timestamps = eventMat.timestamps';
    
    % Note that these licks now are only the ones within the boundaries of the
    % first and last trial (win) - logical.
    eventVar.licks = eventMat.data(:,1:7)';
    
    % Create a single vector with all licking events - logical.
    eventVar.licksAll = eventMat.data(:,1)' | eventMat.data(:,2)' | ...
        eventMat.data(:,3)' | eventMat.data(:,4)' | eventMat.data(:,5)' | ...
        eventMat.data(:,6)' | eventMat.data(:,7)';
    
    % Create a vector with licking events according to port.
    eventVar.licksPorts = eventMat.data(:,1:7)';
    for i = 1:7
        indices = find(eventVar.licks(i,:) == 1);
        eventVar.licksPorts(i,indices) = i;
    end
    
    eventVar.licksPortsAll = zeros(size(eventVar.licksAll));
    for i = 1:7
        indices = find(eventVar.licks(i,:) == 1);
        eventVar.licksPortsAll(indices) = ...
            eventVar.licksPorts(i,indices);
    end
    
    % Create a vector with the id of the lick port only at the index where it
    % first appears anew. 
    eventVar.licksPortsAllSingle = zeros(size(eventVar.licksAll));
    
    licksIdx = find(eventVar.licksPortsAll ~= 0);
    difference = diff(eventVar.licksPortsAll(licksIdx));
    singlesIdx = licksIdx(find(difference ~=0) + 1);
    
    eventVar.licksPortsAllSingle(singlesIdx) = eventVar.licksPortsAll(singlesIdx);
    eventVar.licksPortsAllSingle(1) = 1;
    
    % Create a vector with single licks when they first appear - logical. 
    eventVar.lickEvents = zeros(size(eventVar.licksAll));
    eventVar.lickEvents(singlesIdx) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cellID_licks = find(tuned_licks == 1);
    % cellID_licks = find(tuned_licks == 1 & tuned_relDistStop == 1);
    numTuned = length(cellID_licks);
    meanFiring = zeros(numTuned, win);
    mi_licks1 = zeros(numTuned, 1);
    mi_licks2 = zeros(numTuned, 1);
    mi_licks3 = zeros(numTuned, 1);
    mi_licks4 = zeros(numTuned, 1);
    mi_licks5 = zeros(numTuned, 1);
    
    % Group cells tuned to licks based on which lick provides the most MI.
    best_lick = zeros(numTuned, 1);
    
    j = 0;  % counts the number of cells that are tuned & have valid MI 
    
    for c = 1:numTuned
        j = j + 1;
    
     % (Option 1) Raw firing rate from the recording
        times_spikes = spk_timestamps(spkMat(pyramidal(cellID_licks(c)),:) > 0);
        times_licks = event_timestamps(eventVar.lickEvents == 1);
    
        [ccg,t] = CCG({times_spikes' times_licks'}, [],'binSize',dt, ...
            'duration',2,'norm','rate');
        psth = ccg(1:end-1,2,1)'; 
    %     psth = ccg(:,2,1)'; 
        meanFiring(j, :) = psth;
    %         figure
    %         plot(t,psth)
    
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
                string(pyramidal(cellID_licks(c))), '_fit_k-fold.csv')))
            results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
                    string(pyramidal(cellID_licks(c))), '_fit_k-fold.csv')));
        else
            continue
        end
        % FIXME this might be what is classified as DISCARD NEURON?
        % Ensure that there are values for the MI 
        if ~isnan(table2array(results(4, 11))) & ~isnan(table2array(results(5, 11))) ...
                & ~isnan(table2array(results(6, 11))) & ~isnan(table2array(results(7, 11))) ...
                & ~isnan(table2array(results(8, 11)))
    
            % (Option 2) Raw firing rate from the PGAM (15 values avail -> win=15)
    %         raw_rate = strtrim(results.raw_rate_Hz{11});
    %         meanFiring(j, :) = str2double(regexp(raw_rate, '\d+\.?\d*', 'match'));
    
            % Find MI for each type of lick 
            mi_licks1(j, 1) = table2array(results(4, 11));
            mi_licks2(j, 1) = table2array(results(5, 11));
            mi_licks3(j, 1) = table2array(results(6, 11));
            mi_licks4(j, 1) = table2array(results(7, 11));
            mi_licks5(j, 1) = table2array(results(8, 11));
    
            % Find type of lick with highest MI
            mi_best = max([table2array(results(4, 11)), table2array(results(5, 11)), ...
            table2array(results(6, 11)), table2array(results(7, 11)), ...
            table2array(results(8, 11))]);
    
            rowIndex = find(results.mutual_info == mi_best);
            best_lick(j, 1) = rowIndex - 5;
        else
            continue
        end  
    end
    
    num_tunedMI = 1:j;
    
    % Sort cells based on when max firing occurs
    [~,idx] = max(meanFiring(:,:), [], 2);
    [~,sortidx] = sort(idx);
    meanFiring_sorted = meanFiring(sortidx, :);
    
    figure;
    sp1 = subplot(1,2,1);
    imagesc(zscore(meanFiring_sorted,[],2));
    colorbar
    xline(win/2, 'LineWidth', 1, 'Color', 'r');
    xlabel('Time (2 s)')
    ylabel('Neuron')
    title('Firing around the licking event for all neurons tuned to licks')
    
    sp2 = subplot(1,2,2);
    set(gca, 'YDir','reverse')
    hold on
    scatter(ones(size(num_tunedMI(best_lick==1)))*0.25, num_tunedMI(best_lick==1), 25, [0 0.4470 0.7410],'filled');
    scatter(ones(size(num_tunedMI(best_lick==2)))*0.5, num_tunedMI(best_lick==2), 25, [0.8500 0.3250 0.0980],'filled');
    scatter(ones(size(num_tunedMI(best_lick==3)))*0.75, num_tunedMI(best_lick==3), 25, [0.9290 0.6940 0.1250],'filled');
    scatter(ones(size(num_tunedMI(best_lick==4)))*1, num_tunedMI(best_lick==4), 25, [0.4940 0.1840 0.5560],'filled');
    scatter(ones(size(num_tunedMI(best_lick==5)))*1.25, num_tunedMI(best_lick==5), 25, [0.4660 0.6740 0.1880],'filled');
    hold off
    xticks([0.25, 0.5, 0.75, 1, 1.25])
    xticklabels(["licks1", "licks2", "licks3", "licks4", "licks5"])
    xlim([0,1.5])
    title('Lick type tuning specificity based on max Mutual Information')
    
    saveas(gcf, fullfile(postprocessDir, 'licksTuning_heatmap'), 'fig');
    saveas(gcf, fullfile(postprocessDir, 'licksTuning_heatmap'), 'png');
    close(figure(1));
    
    %% Analyze MI for different lick types
    
    % Plot distributions of MI 
    mi_licks = horzcat(mi_licks1, mi_licks2, mi_licks3, mi_licks4, mi_licks5);
    mi_licks_norm = (mi_licks - min(mi_licks(:))) / (max(mi_licks(:)) - min(mi_licks(:)));
    
    figure;
    set(gcf, 'Position', [100, 100, 800, 600]);
    hold on
    
    map = colormap('hot');
    indices = round(linspace(1, 180, 5));
    colors = map(indices,:);
    
    for i = 1:5
        j = i + 5;
        subplot(2,5,i)
        histogram(mi_licks(:,i), numTuned, 'FaceColor', colors(i,:), ...
            'EdgeColor', colors(i,:));
        title(sprintf('MI licks %d', i));
        xlim([0,max(mi_licks(:,i))])
    
        subplot(2,5,j)
        histogram(mi_licks_norm(:,i), numTuned, 'FaceColor', colors(i,:), ...
            'EdgeColor', colors(i,:));
        title(sprintf('Normalized MI licks %d', i));
        xlim([0,1])
    end
    hold off
    sgtitle('Distribution of MI for cells tuned to licks')
    
    saveas(gcf, fullfile(postprocessDir, 'MI_licks_distribution'), 'fig');
    saveas(gcf, fullfile(postprocessDir, 'MI_licks_distribution'), 'png');
    close(figure(2));
    
    % % Display MI for different licks
    % thres = 0.1; 
    % % mi_licks_dig = zeros(size(mi_licks_norm));
    % % mi_licks_dig(mi_licks_norm > thres) = 1;
    % 
    % figure;
    % [~,idx] = sort(mi_licks_norm(:,1));
    % imagesc(mi_licks_norm(idx,:),[],2);
    % % colorbar('Ticks',[0,1]);
    % % set(gca, 'clim', [0,1])
    % xticks([1,2,3,4,5]);
    % xlabel('Lick type')
    % ylabel('Cell')
    % title(sprintf('Tuning to different lick types with MI thres = %.1f', thres));
    % 
    % saveas(gcf, fullfile(postprocessDir, 'MI_licks_heatmap'), 'fig');
    % saveas(gcf, fullfile(postprocessDir, 'MI_licks_heatmap'), 'png');
    % close(figure(3));
    
    %% Scatter plot of max MI for lick types 1, 4, 5
    colors = zeros(length(cellID_licks), 3);
    maxlick = zeros(length(cellID_licks), 1);
    
    % First normalize and then find the max
    mi_licks1 = normalize(mi_licks1, "range");
    mi_licks4 = normalize(mi_licks4, "range");
    mi_licks5 = normalize(mi_licks5, "range");
    
    % Check that the file exists, otherwise set the max MI value to NaN
    for i = 1:length(cellID_licks)
        if ~exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
                string(cellID_licks(i)), '_fit_k-fold.csv')))
            mi_licks1(i,1) = NaN;
            mi_licks4(i,1) = NaN;
            mi_licks5(i,1) = NaN;
            maxlick(i,1) = NaN;
        else
    %         maxlick(i,1) = max([mi_licks1(i), mi_licks4(i), mi_licks5(i)]);
            maxlick(i,1) = max([mi_licks1(i), mi_licks4(i)]);
        end
    
        % Here we are setting the max MI to be that of lick5 for cells that 
        % have the same MI (usually 0) for more than one lick type e.g. both
        % are max (=1)
        if ((maxlick(i,1) == mi_licks1(i)) && (maxlick(i,1) == mi_licks4(i)))  %| ...
    %             ((maxlick(i,1) == mi_licks1(i)) && (maxlick(i,1) == mi_licks5(i))) | ...
    %             ((maxlick(i,1) == mi_licks4(i)) && (maxlick(i,1) == mi_licks5(i))) | ...
    %             ((maxlick(i,1) == mi_licks1(i)) && (maxlick(i,1) == mi_licks4(i)) && ...
    %             (maxlick(i,1) == mi_licks5(i)))
    %         maxlick(i,1) = mi_licks5(i);
            maxlick(i,1) = mi_licks1(i);
        end
        
        if ~isnan(maxlick(i,1))
            if maxlick(i,1) == mi_licks1(i,1)
                colors(i,:) = [0.6350 0.0780 0.1840];
            elseif maxlick(i,1) == mi_licks4(i,1)
                colors(i,:) = [0 0.4470 0.7410];
            else 
                colors(i,:) = [0.4660 0.6740 0.1880];
            end
        end
    end
    
    figure;
    hold on
    grid on
    % scatter3(mi_licks1(maxlick==mi_licks1), mi_licks4(maxlick==mi_licks1), ...
    %     mi_licks5(maxlick==mi_licks1), [], colors(maxlick==mi_licks1,:), ...
    %     'filled', 'DisplayName', 'Max is licks1');
    % scatter3(mi_licks1(maxlick==mi_licks4), mi_licks4(maxlick==mi_licks4), ...
    %     mi_licks5(maxlick==mi_licks4), [], colors(maxlick==mi_licks4,:), ...
    %     'filled', 'DisplayName', 'Max is licks4');
    % scatter3(mi_licks1(maxlick==mi_licks5), mi_licks4(maxlick==mi_licks5), ...
    %     mi_licks5(maxlick==mi_licks5), [], colors(maxlick==mi_licks5,:), ...
    %     'filled', 'DisplayName', 'Max is licks5');
    scatter(mi_licks1(maxlick==mi_licks1), mi_licks4(maxlick==mi_licks1), [], ...
        colors(maxlick==mi_licks1,:), 'filled', 'DisplayName', 'Max is licks1')
    scatter(mi_licks1(maxlick==mi_licks4), mi_licks4(maxlick==mi_licks4), [], ...
        colors(maxlick==mi_licks4,:), 'filled', 'DisplayName', 'Max is licks4')
    
    xlabel('MI licks1')
    ylabel('MI licks4')
    % zlabel('MI licks5')
    legend
    saveas(gcf, fullfile(postprocessDir, 'licksTuning_scatter2D'), 'fig');
    saveas(gcf, fullfile(postprocessDir, 'licksTuning_scatter2D'), 'png');
    
    %% Calculate proportion of cells tuned to each variable. 
    
    % All cells (pyramidal + interneurons)
    tuned_relDistStop_all = zeros(numNeurons, 1);
    
    for i = 1:numNeurons
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(i), '_fit_k-fold.csv')))
    
            results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
                string(i), '_fit_k-fold.csv')));
        else
            continue
        end
    
        if table2array(results(3, 10)) < 0.001
            tuned_relDistStop_all(i, 1) = 1; 
        end
    end
    
    resultsPGAM_UMAP = {};
    resultsPGAM_UMAP.tuned_relDistStop = tuned_relDistStop_all;
    
    % Only pyramidal cells: this changes the indices. 
    tuned_y = zeros(length(pyramidal),1);
    tuned_ylin = zeros(length(pyramidal),1);
    tuned_relDistStop = zeros(length(pyramidal),1);
    
    % Find the cells tuned to the rest of the variables
    for i = 1:length(pyramidal)
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(pyramidal(i)), '_fit_k-fold.csv')))
            results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
                string(pyramidal(i)), '_fit_k-fold.csv')));
        else
            continue
        end
    
        if table2array(results(1, 10)) < 0.001
            tuned_y(i, 1) = 1; 
        end
        if table2array(results(2, 10)) < 0.001
            tuned_ylin(i, 1) = 1; 
        end
        if table2array(results(3, 10)) < 0.001
            tuned_relDistStop(i, 1) = 1; 
        end
    end
    
    resultsPGAM.tunedY = tuned_y;
    resultsPGAM.tunedYlin = tuned_ylin;
    resultsPGAM.tunedrelDistStop = tuned_relDistStop;
    
    % Find the proportion of all cells tuned to each variable
    yProp = length(find(tuned_y == 1)) / length(pyramidal);
    ylinProp = length(find(tuned_ylin == 1)) / length(pyramidal);
    relDistStopProp = length(find(tuned_relDistStop == 1)) / length(pyramidal);
    licksProp = numTuned / length(pyramidal);
    
    resultsPGAM.yProp = yProp;
    resultsPGAM.ylinProp = ylinProp;
    resultsPGAM.relDistStopProp = relDistStopProp;
    resultsPGAM.licksProp = licksProp;
    
    cells_y_relDistStop = intersect(find(tuned_y == 1), find(tuned_relDistStop == 1));
    cells_licks_relDistStop = intersect(find(tuned_relDistStop == 1), find(tuned_licks == 1));
    cells_y_licks = intersect(find(tuned_licks == 1), find(tuned_y == 1));
    cells_all = find(tuned_y == 1 & tuned_licks == 1 & tuned_relDistStop == 1);
    
    y_relDistStop_prop = length(cells_y_relDistStop) / length(pyramidal);
    licks_relDistStop_prop = length(cells_licks_relDistStop) / length(pyramidal);
    y_licks_prop = length(cells_y_licks) / length(pyramidal);
    y_licks_relDistStop_prop = length(cells_all) / length(pyramidal);
    
    variables = [yProp, ylinProp, relDistStopProp, licksProp, ...
        y_relDistStop_prop, licks_relDistStop_prop, y_licks_prop, ...
        y_licks_relDistStop_prop];
    variable_names = ["y", "ylin", "relDistStop", "licks", "y relDistStop", ...
        "licks relDistStop", "y licks", "y licks relDistStop"];
    
    % Plot the proportions 
    figure;
    hold on
    for v = 1:length(variables)
        bar(v, variables(v), 'EdgeColor', [0.4 0.4 0.4], 'FaceColor', [0.7 0.7 0.7])
    %     scatter(v, variables(v), 'filled', 'DisplayName', sprintf(variable_names(v)))
    end
    xlim([0, length(variables)+1])
    ylim([0,1])
    set(gca,'XTick', 1:length(variables), 'XTickLabel', variable_names);
    % xticklabels(variable_names)
    % legend
    
    saveas(gcf, fullfile(postprocessDir, 'percVariableTuning'), 'fig');
    saveas(gcf, fullfile(postprocessDir, 'percVariableTuning'), 'png');
    
    %% Determine if cells tuned to relDistStop are also tuned to ylin or y 
    % and vice versa. Also look at how the mutual information is distributed.
    % Note that the indices here refer again to the pyramidal cells, not the
    % absolute cell ID.
    
    % % ylin vs relDistStop 
    % % cells_ylin_relDistStop = find(tuned_ylin == 1 | tuned_relDistStop == 1);
    % cells_ylin_relDistStop = intersect(find(tuned_ylin == 1), find(tuned_relDistStop == 1));
    % 
    % % bothPropA = length(cells_ylin_relDistStop) / length(pyramidal);
    % % ylin_bothA = length(cells_ylin_relDistStop) / length(find(tuned_ylin == 1));
    % % relDistStop_bothA = length(cells_ylin_relDistStop) / length(find(tuned_relDistStop == 1));
    % 
    % mi_ylin = zeros(length(cells_ylin_relDistStop), 1);
    % mi_relDistStop = zeros(length(cells_ylin_relDistStop), 1);
    % 
    % for i = 1:length(cells_ylin_relDistStop)
    %     if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %         string(pyramidal(cells_ylin_relDistStop(i))), '_fit_k-fold.csv')))
    % 
    %         results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %             string(pyramidal(cells_ylin_relDistStop(i))), '_fit_k-fold.csv')));
    %         
    %         mi_ylin(i, 1) = table2array(results(3, 11));
    %         mi_relDistStop(i, 1) = table2array(results(4, 11));
    %         
    %     else
    %         continue
    %     end
    % end
    % 
    % mi_relDistStop = normalize(mi_relDistStop, "range");
    % mi_ylin = normalize(mi_ylin, "range");
    % 
    % figure;
    % subplot(1,3,1)
    % scatter(mi_relDistStop, mi_ylin, 'filled')
    % hold on
    % plot([0,max([mi_relDistStop, mi_ylin])], [0,max([mi_relDistStop, mi_ylin])]);
    % xlabel('relDistStop Mutual Information');
    % ylabel('ylin Mutual Information');
    % title('Mutual Information for ylin and relDistStop for cells tuned to both.')
    % hold off
    
    % y vs relDistStop 
    cells_y_relDistStop = pyramidal(tuned_y == 1 & tuned_relDistStop == 1);
    
    mi_y = zeros(length(cells_y_relDistStop), 1);
    mi_relDistStop = zeros(length(cells_y_relDistStop), 1);
    
    tuned_relDistStop_mi = zeros(length(cells_y_relDistStop), 1);
    
    for i = 1:length(cells_y_relDistStop)
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(cells_y_relDistStop(i)), '_fit_k-fold.csv')))
    
            results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(cells_y_relDistStop(i)), '_fit_k-fold.csv')));
    
            mi_y(i, 1) = table2array(results(2, 11));
            mi_relDistStop(i, 1) = table2array(results(4, 11));
    
        else
            continue
        end   
    end
    
    mi_relDistStop = normalize(mi_relDistStop, "range");
    mi_y = normalize(mi_y, "range");
    
    for i = 1:length(cells_y_relDistStop)
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(cells_y_relDistStop(i)), '_fit_k-fold.csv')))
    
            if mi_y(i, 1) < mi_relDistStop(i, 1)
                tuned_relDistStop_mi(i, 1) = 1;
            end
        else
            continue
        end   
    end
    
    resultsPGAM.tuned_relDistStop_mi = cells_y_relDistStop(tuned_relDistStop_mi==1)';
    resultsPGAM.tuned_relDistStop_mi_licks = intersect(resultsPGAM.tuned_relDistStop_mi, ...
        pyramidal(tuned_licks==1));
    resultsPGAM.tuned_relDistStop_mi_nolicks = intersect(resultsPGAM.tuned_relDistStop_mi, ...
        pyramidal(tuned_licks==0));
    
    
    % subplot(1,3,2)
    % scatter(mi_relDistStop, mi_y, 'filled')
    % hold on
    % plot([0,max([mi_relDistStop, mi_y])], [0,max([mi_relDistStop, mi_y])]);
    % xlabel('relDistStop Mutual Information');
    % ylabel('y Mutual Information');
    % title('Mutual Information for y and relDistStop for cells tuned to both.')
    % hold off
    % 
    % % y vs ylin 
    % % cells_y_ylin = find(tuned_y == 1 | tuned_ylin == 1);
    % cells_y_ylin = intersect(find(tuned_y == 1), find(tuned_ylin == 1));
    % 
    % % bothPropC = length(cells_y_ylin) / length(pyramidal);
    % % y_bothC = length(cells_y_ylin) / length(find(tuned_y == 1));
    % % ylin_bothC = length(cells_y_ylin) / length(find(tuned_ylin == 1));
    % 
    % mi_y = zeros(length(cells_y_ylin), 1);
    % mi_ylin = zeros(length(cells_y_ylin), 1);
    % 
    % for i = 1:length(cells_y_ylin)
    %     if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %         string(pyramidal(cells_y_ylin(i))), '_fit_k-fold.csv')))
    % 
    %         results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %             string(pyramidal(cells_y_ylin(i))), '_fit_k-fold.csv')));
    % 
    %         mi_y(i, 1) = table2array(results(2, 11));
    %         mi_ylin(i, 1) = table2array(results(3, 11));
    % 
    %     else 
    %         continue
    %     end
    % end
    % 
    % mi_ylin = normalize(mi_ylin, "range");
    % mi_y = normalize(mi_y, "range");
    % 
    % subplot(1,3,3)
    % scatter(mi_ylin, mi_y, 'filled')
    % hold on
    % plot([0,max([mi_ylin, mi_y])], [0,max([mi_ylin, mi_y])]);
    % xlabel('ylin Mutual Information');
    % ylabel('y Mutual Information');
    % title('Mutual Information for y and ylin for cells tuned to both.')
    % hold off
    % 
    % saveas(gcf, fullfile(postprocessDir, 'y_ylin_relDistStop_MI'), 'fig');
    % saveas(gcf, fullfile(postprocessDir, 'y_ylin_relDistStop_MI'), 'png');
    % 
    % % Examine how MI is distributed for cells tuned to all three variables. 
    % cells_all = find(tuned_y == 1 & tuned_ylin == 1 & tuned_relDistStop == 1);
    % cells_y_ylin_only = find(tuned_y == 1 & tuned_ylin == 1 & tuned_relDistStop == 0);
    % cells_all_id = zeros(length(cells_all), 1);
    % cells_y_ylin_only_id = zeros(length(cells_y_ylin_only), 1);
    % 
    % mi_y_all = zeros(length(cells_all), 1);
    % mi_ylin_all = zeros(length(cells_all), 1);
    % mi_relDistStop_all = zeros(length(cells_all), 1);
    % 
    % mi_y_part = zeros(length(cells_y_ylin_only), 1);
    % mi_ylin_part = zeros(length(cells_y_ylin_only), 1);
    % mi_relDistStop_part = zeros(length(cells_y_ylin_only), 1);
    % 
    % for i = 1:length(cells_all)
    %     if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %         string(pyramidal(cells_all(i))), '_fit_k-fold.csv')))
    % 
    %         results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %             string(pyramidal(cells_all(i))), '_fit_k-fold.csv')));
    %     
    %         mi_y_all(i, 1) = table2array(results(2, 11));
    %         mi_ylin_all(i, 1) = table2array(results(3, 11));
    %         mi_relDistStop_all(i, 1) = table2array(results(4, 11));
    %     
    %         cells_all_id(i, 1) = cells_all(i);
    % 
    %     else
    %         continue
    %     end
    % end
    % 
    % cells_all_id_thres = cells_all_id(mi_relDistStop_all > 0.8);
    % 
    % for i = 1:length(cells_y_ylin_only)
    %     if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %         string(pyramidal(cells_y_ylin_only(i))), '_fit_k-fold.csv')))
    % 
    %         results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
    %             string(pyramidal(cells_y_ylin_only(i))), '_fit_k-fold.csv')));
    %     
    %         mi_y_part(i, 1) = table2array(results(2, 11));
    %         mi_ylin_part(i, 1) = table2array(results(3, 11));
    %         mi_relDistStop_part(i, 1) = table2array(results(4, 11));
    %     
    %         cells_y_ylin_only_id(i, 1) = cells_y_ylin_only(i);
    % 
    %     else
    %         continue
    %     end
    % end
    % 
    % mi_y_all = normalize(mi_y_all, "range");
    % mi_ylin_all = normalize(mi_ylin_all, "range");
    % mi_relDistStop_all = normalize(mi_relDistStop_all, "range");
    % mi_y_part = normalize(mi_y_part, "range");
    % mi_ylin_part = normalize(mi_ylin_part, "range");
    % mi_relDistStop_part = normalize(mi_relDistStop_part, "range");
    % 
    % figure;
    % scatter3(mi_y_all, mi_ylin_all, mi_relDistStop_all, 'filled', 'DisplayName', ...
    %     'Cells tuned to y & ylin & relDistStop');
    % xlabel('MI y');
    % ylabel('MI ylin')
    % zlabel('MI relDistStop')
    % hold on 
    % 
    % scatter3(mi_y_part, mi_ylin_part, mi_relDistStop_part, 'filled', 'DisplayName', ...
    %     'Cells tuned to y & ylin');
    % legend
    % 
    % saveas(gcf, fullfile(postprocessDir, 'y_ylin_relDistStop_scatter'), 'fig');
    % saveas(gcf, fullfile(postprocessDir, 'y_ylin_relDistStop_scatter'), 'png');
    
    %% Find cells tuned to combinations of variables
    tuned_relDistStop_licks = pyramidal(tuned_relDistStop == 1 & tuned_licks == 1)';
    tuned_relDistStop_noLicks = pyramidal(tuned_relDistStop == 1 & tuned_licks == 0)';
    
    resultsPGAM.tuned_relDistStop_licks = pyramidal(tuned_relDistStop == 1 & tuned_licks == 1)';
    resultsPGAM.tuned_relDistStop_noLicks = pyramidal(tuned_relDistStop == 1 & tuned_licks == 0)';
    resultsPGAM.tuned_licks = pyramidal(tuned_licks == 1)';
    resultsPGAM.tuned_y = pyramidal(tuned_y == 1)';
    resultsPGAM.tuned_ylin = pyramidal(tuned_ylin == 1)';
    resultsPGAM.tuned_relDistStop = pyramidal(tuned_relDistStop == 1)';
    
    resultsPGAM.propTuned_relDistStop_licks = length(tuned_relDistStop_licks) / ...
        length(find(tuned_relDistStop==1));
    resultsPGAM.propTuned_relDistStop_noLicks = length(tuned_relDistStop_noLicks) / ...
        length(find(tuned_relDistStop==1));
    
    % Find cells with positive/negative kernel strength to licks
    kernelStrength_licks = zeros(length(pyramidal),1);
    
    for i = 1:length(pyramidal)
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(pyramidal(i)), '_fit_k-fold.csv')))
    
            results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
                string(pyramidal(i)), '_fit_k-fold.csv')));
        else
            continue
        end
    
        if table2array(results(9, 18)) > 0
            kernelStrength_licks(i, 1) = 1; 
        end
    end
    
    resultsPGAM.kernelStrength_licks = pyramidal(kernelStrength_licks==1 & tuned_licks==1)';
    
    %% Plot all their maps together 
    % % rows = ceil(sqrt(length(tuned_relDistStop_licks)));
    % % columns = ceil(sqrt(length(tuned_relDistStop_licks)));
    % 
    % figure('WindowState','maximized','Name','Cells tuned to relDistStop and licks');
    % rows = ceil(length(tuned_relDistStop_licks)/6);
    % columns = ceil(length(tuned_relDistStop_licks)/4);
    % 
    % for i = 1:length(tuned_relDistStop_licks)
    %     cell = tuned_relDistStop_licks(i);
    %     imageName = fullfile(mapDir, strcat('cell_', string(cell), '.png'));
    %     image = imread(imageName);
    % 
    %     subaxis(rows, columns, i, 'SpacingVert',0.01,'SpacingHoriz',0.01, ...
    %         'MarginLeft',0.01,'MarginRight',0.01,'MarginTop',0.01,'MarginBottom',0.01);
    %     imshow(image); 
    %     title(strcat('cell ', string(cell)));
    % end
    % 
    % % sgtitle('Cells tuned to relDistStop and licks');
    % saveas(gcf, fullfile(postprocessDir, 'maps_tuned_relDistStop_licks'), 'fig');
    % 
    % figure('WindowState','maximized','Name','Cells tuned to relDistStop but not to licks');
    % rows = ceil(length(tuned_relDistStop_noLicks)/6);
    % columns = ceil(length(tuned_relDistStop_noLicks)/5);
    % 
    % for i = 1:length(tuned_relDistStop_noLicks)
    %     cell = tuned_relDistStop_noLicks(i);
    %     imageName = fullfile(mapDir, strcat('cell_', string(cell), '.png'));
    %     image = imread(imageName);
    % 
    %     subaxis(rows, columns, i, 'SpacingVert',0.01,'SpacingHoriz',0.01, ...
    %         'MarginLeft',0.01,'MarginRight',0.01,'MarginTop',0.01,'MarginBottom',0.01);
    %     imshow(image);
    %     title(strcat('cell ', string(cell)));
    % end
    % 
    % % sgtitle('Cells tuned to relDistStop but not to licks');
    % saveas(gcf, fullfile(postprocessDir, 'maps_tuned_relDistStop_noLicks'), 'fig');
    
    %% Save results
    if exist(fullfile(postprocessDir, strcat(extractAfter(dataDir, 'Videos\'), ...
        '.postprocessPGAM.mat')))
        delete(fullfile(postprocessDir, strcat(extractAfter(dataDir, 'Videos\'), ...
        '.postprocessPGAM.mat')));
    end
    save(fullfile(postprocessDir, strcat(extractAfter(dataDir, 'Videos\'), ...
        '.postprocessPGAM.mat')), 'resultsPGAM'); 
    
    if exist(fullfile(postprocessDir, strcat(extractAfter(dataDir, 'Videos\'), ...
        '.resultsPGAM_UMAP.mat')))
        delete(fullfile(postprocessDir, strcat(extractAfter(dataDir, 'Videos\'), ...
        '.resultsPGAM_UMAP.mat')));
    end
    save(fullfile(postprocessDir, strcat(extractAfter(dataDir, 'Videos\'), ...
        '.resultsPGAM_UMAP.mat')), 'resultsPGAM_UMAP'); 
end
end