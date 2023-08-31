load('IZ50_230621_sess9.firingMapsTrial.cellinfo.mat')

load('IZ50_230621_sess9.directionStats.cellinfo.mat')
placeCells = unique([directionStats.placeCells.right, directionStats.placeCells.left]);
datamat_mean{1} = [];
datamat_mean{2} = [];

data{1} = [];
data{2} = [];
for pf = 1:length(placeCellsAll)
    datamat{1} = [];
    datamat{2} = [];
    
    for kk = 1:length(firingMaps.right.rateMaps{placeCellsAll(pf)})
    datamat{1} = [datamat{1}; firingMaps.right.rateMaps{placeCellsAll(pf)}{kk}];
    datamat{2} = [datamat{2}; firingMaps.left.rateMaps{placeCellsAll(pf)}{kk}];
    end
    datamat_mean{1} = [datamat_mean{1}; mean(datamat{1})];
    datamat_mean{2} = [datamat_mean{2}; mean(datamat{2})];
    
    data{1} = [data{1}; firingMaps.rateMaps{placeCellsAll(pf)}{1}];
    data{2} = [data{2}; firingMaps.rateMaps{placeCellsAll(pf)}{2}];
end
figure
subplot(2,2,1)
imagesc(datamat_mean{1})
subplot(2,2,2)
imagesc(datamat_mean{2})
subplot(2,2,3)
imagesc(data{1})
subplot(2,2,4)
imagesc(data{2})

%% Misc functions
%% Plot correlation across trials for each cell 
% if ~isfolder('FiringMap/CellTrialFR')
%     mkdir('FiringMap/CellTrialFR');
% end 
% 
% for ll = 1:length(labels)
%     for pf = 1:length(placeCells{ll})
%         cells = placeCells{ll};
%         figure
%         set(gcf,'Renderer','painters')
%         set(gcf,'Position',[2200 200 1185 712])
%         datamat_lin{1} = [];
%         datamat_lin{2} = [];
%         
%         for kk = 1:numtrials
%             datamat_lin{ll} = [datamat_lin{ll}, firingMaps.(labels{ll}).rateMaps{cells(pf)}{kk}];
%         end
%         
%         [corr,lags] = xcorr(datamat_lin{ll});
%         
%         [maxima_values, maxima_locs] = findpeaks(corr,'MinPeakDistance',55); 
%         
%         % Plot the autocorrelation function and the local maxima
%         figure;
%         plot(lags, corr, 'k-');
%         hold on;
%         plot(maxima_locs+lags(1), maxima_values, 'b-', 'Marker', '.')
%         
%         % Fit an exponential model
%         modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1)) + b(3);  
%         beta0 = [0, 0, 0]; 
%         
%         % Create the model and determine the coefficients
% %         tbl = table(maxima_locs(length(maxima_locs)/2:end)'+lags(1), maxima_values(length(maxima_locs)/2:end)');
%         tbl = table(lags(length(lags)/2:end)', yy2);
%         mdl = fitnlm(tbl, modelfun, beta0);
%         coefficients = mdl.Coefficients.Estimate;
% 
%         % Create smoothed/regressed data using the model
% %         yFitted = coefficients(1) * exp(-coefficients(2)*(maxima_locs(length(maxima_locs)/2:end)'+lags(1))) + coefficients(3);
%         yFitted = coefficients(1) * exp(-coefficients(2)*(lags(length(lags)/2:end))) + coefficients(3);
%         
%         plot(lags(length(lags)/2:end), yFitted, 'r-', 'LineWidth', 2);
%         
%         % (Optional) Interpolate
%         yint = interp1(maxima_locs(length(maxima_locs)/2:end)+lags(1), ...
%             maxima_values(length(maxima_locs)/2:end), ... 
%             lags(length(lags)/2:end));
%         yy2 = smooth(lags(length(lags)/2:end),yint,0.05,'lowess');
%         figure
%         plot(lags(length(lags)/2:end), yy2)
%         % Determine the steepness of the fit curve
%         [fy] = gradient(yint);
% 
%         saveas(gcf,['FiringMap/CellTrialFR',filesep,'cell_' num2str(cells(pf)) '_' labels{ll} '.png'],'png');
%         saveas(gcf,['FiringMap/CellTrialFR',filesep,'cell_' num2str(cells(pf)) '_' labels{ll} '.fig'],'fig');
%         close all;
%      end
% end

%% Plot the cells with fields in both directions
% subplots = 2*length(cells2fields);
% 
% if subplots <= 16
%     figure;
%     set(gcf,'Position',[300 100 1300 700]);
% 
%     n = 1;
%     for d = 1:directions
%         for unit = 1:length(cells2fields)
%             subplot(directions, length(cells2fields), n);
%             plot(firingMaps.rateMaps{cells2fields(unit)}{d},'k');
%             if sum(firingMaps.rateMaps{cells2fields(unit)}{d}) > 0
%                 hold on
%                 for ii = 1:size(placeFieldStats.mapStats{cells2fields(unit)}{d}.field, 2)
%                     plot(find(placeFieldStats.mapStats{cells2fields(unit)}{d}.field(:,ii)), ...
%                         firingMaps.rateMaps{cells2fields(unit)}{d}(placeFieldStats.mapStats{cells2fields(unit)}{d}.field(:,ii)==1),'linewidth',2)
%                     plot([1 1]*placeFieldStats.mapStats{cells2fields(unit)}{d}.x(ii), ...
%                         [0 firingMaps.rateMaps{cells2fields(unit)}{d}(placeFieldStats.mapStats{cells2fields(unit)}{d}.x(ii)==1)],'--k')
%                 end
%             end
%             set(gca,'XTickLabel',[]);
%             n = n + 1;
%         end
%         
%         if ~exist([basepath filesep 'newPCs'])
%             mkdir([basepath filesep 'newPCs'])
%         end
%     
%         annotation('textbox', [0.5, 0.95, 0.8, 0.05], 'String', labels{1}, ...
%             'FontSize', 16, 'FontWeight', 'bold', 'LineStyle', 'none');
%         annotation('textbox', [0.5, 0.48, 0.8, 0.05], 'String', labels{2}, ...
%             'FontSize', 16, 'FontWeight', 'bold', 'LineStyle', 'none'); 
%     end
%     saveas(gcf,[basepath,filesep,'newPCs',filesep ,'map_2dirs' '.fig'],'fig');
%     saveas(gcf,[basepath,filesep,'newPCs',filesep ,'map_2dirs' '.png'],'png');
%     close(figure(1))
% 
% else % make two plots 
%     for d = 1:directions
%         figure;
%         set(gcf,'Position',[300 100 1300 700]);
%         for unit = 1:length(cells2fields)
%             subplot(3, ceil(length(cells2fields)/3), unit);
%             plot(firingMaps.rateMaps{cells2fields(unit)}{d},'k');
%             if sum(firingMaps.rateMaps{cells2fields(unit)}{d}) > 0
%                 hold on
%                 for ii = 1:size(placeFieldStats.mapStats{cells2fields(unit)}{d}.field, 2)
%                     plot(find(placeFieldStats.mapStats{cells2fields(unit)}{d}.field(:,ii)), ...
%                         firingMaps.rateMaps{cells2fields(unit)}{d}(placeFieldStats.mapStats{cells2fields(unit)}{d}.field(:,ii)==1),'linewidth',2)
%                     plot([1 1]*placeFieldStats.mapStats{cells2fields(unit)}{d}.x(ii), ...
%                         [0 firingMaps.rateMaps{cells2fields(unit)}{d}(placeFieldStats.mapStats{cells2fields(unit)}{d}.x(ii)==1)],'--k')
%                 end
%             end
%             set(gca,'XTickLabel',[]);
%         end
%         
%         if ~exist([basepath filesep 'newPCs'])
%             mkdir([basepath filesep 'newPCs'])
%         end
%         sgtitle(labels{d})
%         
%         saveas(gcf,[basepath,filesep,'newPCs',filesep ,'map_2dirs_' labels{d} '.fig'],'fig');
%         saveas(gcf,[basepath,filesep,'newPCs',filesep ,'map_2dirs_' labels{d} '.png'],'png');
%     end    
% end
%