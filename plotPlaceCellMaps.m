%% (2) plotPlaceCellMaps (place maps for the cells with fields)
function [datamat_cells_norm] = plotPlaceCellMaps(directionStats, thresFR, firingMaps, tracking, behavTrials, spikes, labels, cond)

if ~isfolder('FiringMapAvg/PlaceCells')
    mkdir('FiringMapAvg/PlaceCells');
end

% Apply firing rate threshold
% dt = mean(diff(tracking.timestamps));
% win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
% spkData = bz_SpktToSpkmat(spikes,'dt',dt,'win',win);
% 
% spkMat = spkData.data';
% timestamps = spkData.timestamps';
% 
% fr = zeros(1, length(spkMat(:,1)));
% for i = 1:length(spkMat(:,1))
%     fr(i) = mean(spkMat(i,:)) * 1/dt;
% end
% LowFiringCells = find(fr < thresFR);

placeCellsAll = unique([directionStats.placeCells.(labels{1}), directionStats.placeCells.(labels{2})]);

datamat_cells{1} = [];
datamat_cells{2} = [];

for pf = 1:length(placeCellsAll)   
%     if ismember(placeCellsAll(pf), LowFiringCells)
%         continue
%     end
    for ll = 1:length(labels)
        datamat_cells{ll} = [datamat_cells{ll}; firingMaps.rateMaps{placeCellsAll(pf)}{ll}];
    end
end

% Normalize firing rate among two directions
datamat_cells_normAll = zeros(size(datamat_cells{1},1), 2*size(datamat_cells{1},2));
for pf = 1:size(datamat_cells{1},1)
    datamat_cells_normAll(pf,:) = normalize([datamat_cells{1}(pf,:), datamat_cells{2}(pf,:)], 'range');
    norm_data = datamat_cells_normAll(pf,:);
    datamat_cells_norm{1}(pf,:) = norm_data(1:length(datamat_cells{1}(pf,:)));
    datamat_cells_norm{2}(pf,:) = norm_data(length(datamat_cells{1}(pf,:))+1:end);
end

figure;
for ll = 1:length(labels)
    [~,idx] = max(datamat_cells_norm{1}, [], 2);
    [~,sortidx] = sort(idx);
    datamat_cellsSorted{ll} = datamat_cells_norm{ll}(sortidx, :);

    subplot(2,2,ll)
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[2200 200 1185 712])
%     imagesc(zscore(datamat_cellsSorted{ll},[],2));
    imagesc(datamat_cellsSorted{ll})
    colorbar
    xlabel('Position bin')
    ylabel('Neuron')
    title(labels{ll})
end
for ll = 1:length(labels)
    [~,idx] = max(datamat_cells_norm{2}, [], 2);
    [~,sortidx] = sort(idx);
    datamat_cellsSorted{ll} = datamat_cells_norm{ll}(sortidx, :);

    subplot(2,2,ll+2)
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[2200 200 1185 712])
%     imagesc(zscore(datamat_cellsSorted{ll},[],2));
    imagesc(datamat_cellsSorted{ll})
    colorbar
    xlabel('Position bin')
    ylabel('Neuron')
    title(labels{ll})
end

if cond ~= 0
    saveas(gcf,['FiringMapAvg/PlaceCells',filesep,'placeCells_' cond '.png'],'png');
    saveas(gcf,['FiringMapAvg/PlaceCells',filesep,'placeCells_' cond 'fig'],'fig');
else
    saveas(gcf,['FiringMapAvg/PlaceCells',filesep,'placeCells.png'],'png');
    saveas(gcf,['FiringMapAvg/PlaceCells',filesep,'placeCells.fig'],'fig');
end
close(figure(1))

end