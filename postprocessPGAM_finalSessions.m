% This function calculates the percentage of cells for the chosen session
% that are tuned to each of the variables that we have chosen to test with
% the PGAM. This is done across all sessions and animals. 
% The function also identifies the cells that are tuned to licks and plots
% their PSTH. 

function [resultsPGAM] = postprocessPGAM_finalSessions(varargin)

close all 

p = inputParser;
addParameter(p,'basepath','Z:\Homes\zutshi01\Recordings\Auditory_Task',@isstr);
addParameter(p,'fits','I:\April 2023 GAM fits\results_IZ\new_vars2',@isstr);
addParameter(p,'data','I:\Videos',@isstr);
addParameter(p,'postprocess','I:\April 2023 GAM fits\postprocess',@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fits = p.Results.fits;
data = p.Results.data;
postprocess = p.Results.postprocess;

animals = dir(postprocess);
animals = animals([animals.isdir]);
animals = animals(~ismember({animals(:).name},{'.','..'}));

names = [];
sessions = {};
for a = 1:length(animals)
    sessions{a} = dir(fullfile(postprocess, animals(a).name));
    sessions{a} = sessions{a}([sessions{a}.isdir]);
    sessions{a} = sessions{a}(~ismember({sessions{a}(:).name},{'.','..'}));

    names = [names, string(animals(a).name)];
end

numSessions = [];
for a = 1:length(animals)
    numSessionsAnimal = length(sessions{a});
    numSessions = [numSessions, numSessionsAnimal];
end

%% Average the proportions of cells tuned across sessions and animals 

% Save the individual values from all sessions
yProp_all = cell(1, length(animals));
ylinProp_all = cell(1, length(animals));
relDistStopProp_all = cell(1, length(animals));
positiveLicksProp_all = cell(1, length(animals));
yFor_relDistStop_prop_all = cell(1, length(animals));
licksP_relDistStop_prop_all = cell(1, length(animals));
y_licksP_prop_all = cell(1, length(animals));
yFor_licksP_relDistStop_prop_all = cell(1, length(animals));

yFor_mi_prop_all = cell(1, length(animals));
relDistStop_mi_prop_all = cell(1, length(animals));
relDistStop_licksP_mi_prop_all = cell(1, length(animals));

for a = 1:length(animals)
    cd(fullfile(postprocess, animals(a).name))
    
    for s = 1:numSessions(a)
        session = sessions{a}(s).name;
        cd(fullfile(postprocess, animals(a).name, session))
        file = dir([pwd filesep '*.postprocessPGAM.mat']);
        load(file.name);

        yProp_all{a} = [yProp_all{a}, resultsPGAM.proportions.yProp];
        ylinProp_all{a} = [ylinProp_all{a}, resultsPGAM.proportions.ylinProp];
        relDistStopProp_all{a} = [relDistStopProp_all{a}, resultsPGAM.proportions.relDistStopProp];
        positiveLicksProp_all{a} = [positiveLicksProp_all{a}, resultsPGAM.proportions.positiveLicksProp];
        yFor_relDistStop_prop_all{a} = [yFor_relDistStop_prop_all{a}, resultsPGAM.proportions.yFor_relDistStop_prop];
        licksP_relDistStop_prop_all{a} = [licksP_relDistStop_prop_all{a}, resultsPGAM.proportions.licksP_relDistStop_prop];
        y_licksP_prop_all{a} = [y_licksP_prop_all{a}, resultsPGAM.proportions.y_licksP_prop];
        yFor_licksP_relDistStop_prop_all{a} = [yFor_licksP_relDistStop_prop_all{a}, resultsPGAM.proportions.yFor_licksp_relDistStop_prop];
    
        yFor_mi_prop_all{a} = [yFor_mi_prop_all{a}, resultsPGAM.MI.yFor_mi_prop];
        relDistStop_mi_prop_all{a} = [relDistStop_mi_prop_all{a}, resultsPGAM.MI.relDistStop_mi_prop];
        % Note: for the variable below, positive kernels are considered
        % although not specified in the struct field name
        relDistStop_licksP_mi_prop_all{a} = [relDistStop_licksP_mi_prop_all{a}, resultsPGAM.MI.relDistStop_licks_mi_prop];
    end
end

% Calculate the mean
yProp_meanAnimal = zeros(length(animals),1);
ylinProp_meanAnimal = zeros(length(animals),1);
relDistStopProp_meanAnimal = zeros(length(animals),1);
positiveLicksProp_meanAnimal = zeros(length(animals),1);
yFor_relDistStop_prop_meanAnimal = zeros(length(animals),1);
licksP_relDistStop_prop_meanAnimal = zeros(length(animals),1);
y_licksP_prop_meanAnimal = zeros(length(animals),1);
yFor_licksP_relDistStop_prop_meanAnimal = zeros(length(animals),1);

yFor_mi_prop_meanAnimal = zeros(length(animals),1);
relDistStop_mi_prop_meanAnimal = zeros(length(animals),1);
relDistStop_licksP_mi_prop_meanAnimal = zeros(length(animals),1);

for a = 1:length(animals)
    yProp_meanAnimal(a) = mean(yProp_all{a});
    ylinProp_meanAnimal(a) = mean(ylinProp_all{a});
    relDistStopProp_meanAnimal(a) = mean(relDistStopProp_all{a});
    positiveLicksProp_meanAnimal(a) = mean(positiveLicksProp_all{a});
    yFor_relDistStop_prop_meanAnimal(a) = mean(yFor_relDistStop_prop_all{a});
    licksP_relDistStop_prop_meanAnimal(a) = mean(licksP_relDistStop_prop_all{a});
    y_licksP_prop_meanAnimal(a) = mean(y_licksP_prop_all{a});
    yFor_licksP_relDistStop_prop_meanAnimal(a) = mean(yFor_licksP_relDistStop_prop_all{a});

    yFor_mi_prop_meanAnimal(a) = mean(yFor_mi_prop_all{a});
    relDistStop_mi_prop_meanAnimal(a) = mean(relDistStop_mi_prop_all{a});
    relDistStop_licksP_mi_prop_meanAnimal(a) = mean(relDistStop_licksP_mi_prop_all{a});
end


%% Plot the proportions

% Plot A: tuning to each variables
figure('Position', [100,96,1501,900]);
hold on

% Scatter plots
colors = [...
     1.0000    0.7333    0.6196; ... 
    0.8510    0.3255    0.0980; ... 
    1.0000    0.8902    0.6314; ... 
    0.9294    0.6941    0.1255; ... 
    0.6667    0.5098    0.7020; ... 
    0.4941    0.1843    0.5569; ... 
    0.6980    0.8196    0.5412; ... 
    0.4667    0.6745    0.1882; ... 
    0.5216    0.6824    0.7882; ... 
    0    0.4471    0.7412; ...
    0.5    0.5    0.5;
    0   0    0];

c = 1;
for a = 1:length(animals)
    x = 0.25;
    variables = [yProp_all{a}; ylinProp_all{a}; relDistStopProp_all{a}; positiveLicksProp_all{a}; ...
        yFor_relDistStop_prop_all{a}; licksP_relDistStop_prop_all{a}; y_licksP_prop_all{a}; ...
        yFor_licksP_relDistStop_prop_all{a}];

    variables_meanAnimal = [yProp_meanAnimal(a), ylinProp_meanAnimal(a), ...
        relDistStopProp_meanAnimal(a), positiveLicksProp_meanAnimal(a), ...
        yFor_relDistStop_prop_meanAnimal(a), licksP_relDistStop_prop_meanAnimal(a), ...
        y_licksP_prop_meanAnimal(a), yFor_licksP_relDistStop_prop_meanAnimal(a)];
    
    color = colors(c, :);
    colorL = colors(c+1, :);
    for v = 1:size(variables, 1) 
        s(a) = scatter(ones(size(variables(v,:)))*x, variables(v,:), 14, color, 'filled');
        xvals = [x-0.03, x+0.03];
        yvals = repmat(variables_meanAnimal(v),1,2);
        pp(a) = plot(xvals, yvals, "-", "Color", colorL, "LineWidth", 2);
        x = x + 0.25;   
    end
    legend(pp, names, 'AutoUpdate', 'off', 'FontSize', 16, 'EdgeColor', 'none')
    c = c + 2;
end

% Boxplots
variables_all = [cell2mat(yProp_all); cell2mat(ylinProp_all); ...
    cell2mat(relDistStopProp_all); cell2mat(positiveLicksProp_all); ...
        cell2mat(yFor_relDistStop_prop_all); cell2mat(licksP_relDistStop_prop_all); 
        cell2mat(y_licksP_prop_all); cell2mat(yFor_licksP_relDistStop_prop_all)];

positions = [0.35, 0.6, 0.85, 1.1, 1.35, 1.6, 1.85, 2.1];
for v = 1:size(variables_all, 1)
    b = boxplot(variables_all(v,:), 'Positions', positions(v), ...
        'Symbol', '', 'Widths', 0.09);    
end

h = findobj(gca,'Tag','Box');
h1=findobj(gca,'Tag','Upper Whisker');
h2=findobj(gca,'Tag','Lower Whisker');
h3=findobj(gca,'Tag','Upper Adjacent Value');
h4=findobj(gca,'Tag','Lower Adjacent Value');
h5=findobj(gca,'Tag','Median');  

for j = 1:length(h)
    set(h(j),'color',[0.7 0.7 0.7])
    set(h1(j),'lineStyle','-','LineWidth',2,'color',[0.7 0.7 0.7]);
    set(h2(j),'lineStyle','-','LineWidth',2,'color',[0.7 0.7 0.7]);
    set(h3(j),'lineStyle','none'); 
    set(h4(j),'lineStyle','none');
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.7 0.7 0.7],'EdgeColor','none');
end

for j = 1:length(h5)
    set(h5(j), 'color', [1 1 1], 'lineWidth', 2);
end
set(gca,'Children',flipud(get(gca,'Children')));

xticks([0.3, 0.55, 0.8, 1.05, 1.3, 1.55, 1.8, 2.05]);
xticklabels(["Position (tone) (y)", "Position (no-tone)", "distance (d)", "licks", "Position (tone) & distance", ...
        "licks & distance", "Position (tone) & licks", "Position (tone) & licks relDistStop"]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel', a,'FontSize', 16);
set(gca, 'TickLength', [0 0]);
xlim([0.2, 2.2])
ylim([0,1])
yticks([0,0.5,1])
ylabel('Fraction of cells')
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0]);  
set(gca, 'LineWidth', 2);
box off

% Save plots
saveas(gcf, fullfile(postprocess, 'percVariableTuning_sessions'), 'fig');
saveas(gcf, fullfile(postprocess, 'percVariableTuning_sessions'), 'png');

close(figure(1))


% Flot B: tuning to both relDistStop and y (forward direction)
figure('Position', [2020,96,670,670]);
hold on

% Scatter plots
c = 1;
for a = 1:length(animals)
    x = 0.25;
    variables = [yFor_mi_prop_all{a}; relDistStop_mi_prop_all{a}];

    variables_meanAnimal = [yFor_mi_prop_meanAnimal(a), relDistStop_mi_prop_meanAnimal(a)];
    
    color = colors(c, :);
    colorL = colors(c+1, :);
    for v = 1:size(variables, 1) 
        f(a) = scatter(ones(size(variables(v,:)))*x, variables(v,:), 14, color, 'filled');
        xvals = [x-0.03, x+0.03];
        yvals = repmat(variables_meanAnimal(v),1,2);
        ll(a) = plot(xvals, yvals, "-", "Color", colorL, "LineWidth", 2);
        x = x + 0.25;   
    end
    legend(ll, names, 'AutoUpdate', 'off', 'FontSize', 16, 'EdgeColor', 'none')
    c = c + 2;
end

% Boxplots
variables_all = [cell2mat(yFor_mi_prop_all); cell2mat(relDistStop_mi_prop_all)];

positions = [0.35, 0.6];
for v = 1:size(variables_all, 1)
    b = boxplot(variables_all(v,:), 'Positions', positions(v), ...
        'Symbol', '', 'Widths', 0.09);
end

h = findobj(gca,'Tag','Box');
h1=findobj(gca,'Tag','Upper Whisker');
h2=findobj(gca,'Tag','Lower Whisker');
h3=findobj(gca,'Tag','Upper Adjacent Value');
h4=findobj(gca,'Tag','Lower Adjacent Value');
h5=findobj(gca,'Tag','Median');  

for j = 1:length(h)
    set(h(j),'color',[0.7 0.7 0.7])
    set(h1(j),'lineStyle','-','LineWidth',2,'color',[0.7 0.7 0.7]);
    set(h2(j),'lineStyle','-','LineWidth',2,'color',[0.7 0.7 0.7]);
    set(h3(j),'lineStyle','none'); 
    set(h4(j),'lineStyle','none');
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.7 0.7 0.7],'EdgeColor','none');
end

for j = 1:length(h5)
    set(h5(j), 'color', [1 1 1], 'lineWidth', 2);
end
set(gca,'Children',flipud(get(gca,'Children')));

xticks([0.3, 0.55, 0.8])
xticklabels(["y > d", "d > y"])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel', a,'FontSize', 16);
set(gca, 'TickLength', [0 0]);
xlim([0.2, 0.7])
ylim([0,1])
yticks([0,0.5,1])
ylabel('Fraction of cells')
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0]);  
set(gca, 'LineWidth', 2);
box off

saveas(gcf, fullfile(postprocess, 'percVariableTuning_prefSessions'), 'fig');
saveas(gcf, fullfile(postprocess, 'percVariableTuning_prefSessions'), 'png');

close(figure(1));
   
end
