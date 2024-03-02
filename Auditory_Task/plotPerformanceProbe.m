function [outputArg1,outputArg2] = plotPerformanceProbe(varargin)

%% Plot changes in performance in probe vs non-probe trials 

p = inputParser;
addParameter(p,'rootDir',pwd,@isstr); % Folder with all sessions
parse(p,varargin{:});
rootDir = p.Results.rootDir;

contents = dir(rootDir);
probeGroup = contents([contents.isdir]);
probeGroup = probeGroup(~ismember({probeGroup(:).name},{'.','..'}));
probeGroup = probeGroup(contains(string({probeGroup(:).name}),'_probe_'));

cutoffFreq = zeros(length(probeGroup), 1);
for i = 1:length(probeGroup)
    cutoffFreq(i) = str2double(extractBetween(string({probeGroup(i).name}),'probe_','kHz')); % kHz
end

% size = #animals x #frequencies x #probeGroup/freq/animal
accuracy_noProbe = zeros(length(probeGroup), 3);
accuracy_probe = zeros(length(probeGroup), 3);

% Load files and save performance data 
for i = 1:length(probeGroup)
    cd(string({probeGroup(i).name})) 
    sessions = dir(pwd);
    sessions = sessions([sessions.isdir]);
    sessions = sessions(~ismember({sessions(:).name},{'.','..'}));

    for j = 1:length(sessions)
        cd(string({sessions(j).name}))
        cur_sess = pwd;
        
        if ~isempty(dir([cur_sess filesep '*.probeTrials.mat']))
            file = dir([cur_sess filesep '*.probeTrials.mat']);
            load(fullfile(cur_sess, file.name));
        end

        accuracy_noProbe(i, j) = probeTrials.accuracy_noProbe;
        accuracy_probe(i, j) = probeTrials.accuracy_Probe;
        cd ..
    end
    cd ..    
end

% Stack the results from all animals together for each frequency
accuracy_probe_new = [accuracy_probe(2,:), accuracy_probe(5,:); ...
    accuracy_probe(3,:), accuracy_probe(6,:); ...
    accuracy_probe(1,:), accuracy_probe(4,:);];

accuracy_noProbe_new = [accuracy_noProbe(2,:), accuracy_noProbe(5,:); ...
    accuracy_noProbe(3,:), accuracy_noProbe(6,:); ...
    accuracy_noProbe(1,:), accuracy_noProbe(4,:);];

% Statistical test (paired t-test)
for i = 1:size(accuracy_probe_new,1)
    [h, p, ci, stats] = ttest2(accuracy_probe_new(i,:), accuracy_noProbe_new(i,:));
    fprintf('p-value: %.6f\n', p);
end

% Plot the performance
cutoffFreq = unique(cutoffFreq)';

figure 

xvals1 = [1, 2, 3]';
xvals2 = xvals1 + 0.5;
p1 = scatter(repmat(xvals1, 1, size(accuracy_noProbe_new, 2)), accuracy_noProbe_new, ...
    'black', 'filled', 'DisplayName', 'Normal trials');
hold on
p2 = scatter(repmat(xvals2, 1, size(accuracy_probe_new, 2)), accuracy_probe_new, ...
    'blue', 'filled', 'DisplayName', 'Probe trials');
xlim([0.5, 4])
xticks([1.25, 2.25, 3.25])
xticklabels(["4 kHz", "8 kHz", "16 kHz"])
xlabel('Cutoff frequency')
ylabel('Accuracy')
legend([p1(1); p2(1)], 'Location','southeast')

for i = 1:length(cutoffFreq)
    for j = 1:size(accuracy_noProbe_new, 2)
        plot([xvals1(i), xvals2(i)], ...
            [accuracy_noProbe_new(i,j), accuracy_probe_new(i,j)], '-', ...
            'color', [0.5 0.5 0.5],'HandleVisibility','off')
    end
end

saveas(gcf, fullfile(rootDir, strcat('performance_probe')), 'fig');
saveas(gcf, fullfile(rootDir, strcat('performance_probe')), 'png');

%% Plot the distance to target for all animals 
figure('Units', 'centimeters', 'PaperUnits', 'centimeters', 'PaperSize', ...
    [21, 29.7], 'PaperPosition', [1, 1, 21, 29.7]);
hold on
p = 1;
box_pos = [0.25, 0.9, 0.5, 0.05];

for i = 1:length(cutoffFreq)
    freq = strcat(string(cutoffFreq(i)), 'kHz');
    animals = probeGroup(contains(string({probeGroup(:).name}),string(freq)));
    
    portsDiff_probeBlockA_all = [];
    portsDiff_probeBlockB_all = [];
    portsDiff_noProbeBlockA_all = [];
    portsDiff_noProbeBlockB_all = [];

    idxProbe_blockA_all = [];
    idxProbe_blockB_all = [];
    idxNoProbe_blockA_all = [];
    idxNoProbe_blockB_all = [];

    for k = 1:length(animals)
        cd(string({animals(k).name}))
        sessions = dir(pwd);
        sessions = sessions([sessions.isdir]);
        sessions = sessions(~ismember({sessions(:).name},{'.','..'}));
    
        for j = 1:length(sessions)
            cd(string({sessions(j).name}))
            cur_sess = pwd;
        
            if ~isempty(dir([cur_sess filesep '*.TrialBehavior.Events.mat']))
                file = dir([cur_sess filesep '*.TrialBehavior.Events.mat']);
                load(fullfile(cur_sess, file.name));
            end
    
            idxProbe = find(behavTrials.probeTrial==1);
            idxNoProbe = find(behavTrials.probeTrial==0);
            toneGain = behavTrials.toneGain + 1;
            lickLoc = behavTrials.lickLoc + 1; 
    
            % Distance to target for probe trials. Block A = 1-3, Block B = 4-6
            idxProbe_blockA = find(behavTrials.toneTrial==0 ...
            & behavTrials.linTrial==0 & behavTrials.probeTrial==1 & ...
            (toneGain == 1 | toneGain == 2 | toneGain == 3));
        
            idxProbe_blockB = find(behavTrials.toneTrial==0 ...
                & behavTrials.linTrial==0 & behavTrials.probeTrial==1 & ...
                (toneGain == 4 | toneGain == 5 | toneGain == 6));
           
            portsDiff_probeBlockA = toneGain(idxProbe_blockA(:)) - ...
                lickLoc(idxProbe_blockA(:));
            portsDiff_probeBlockB = toneGain(idxProbe_blockB(:)) - ...
                    lickLoc(idxProbe_blockB(:));
            
             % Distance to target for non-probe trials. Block A = 1-3, Block B = 4-6
            idxNoProbe_blockA = find(behavTrials.toneTrial==0 ...
                & behavTrials.linTrial==0 & behavTrials.probeTrial==0 & ...
                (toneGain == 1 | toneGain == 2 | toneGain == 3));
            
            idxNoProbe_blockB = find(behavTrials.toneTrial==0 ...
                & behavTrials.linTrial==0 & behavTrials.probeTrial==0 & ...
                (toneGain == 4 | toneGain == 5 | toneGain == 6));
      
            portsDiff_noProbeBlockA = toneGain(idxNoProbe_blockA(:)) - ...
                lickLoc(idxNoProbe_blockA(:));
            portsDiff_noProbeBlockB = toneGain(idxNoProbe_blockB(:)) - ...
                lickLoc(idxNoProbe_blockB(:));
    
            portsDiff_probeBlockA_all = [portsDiff_probeBlockA_all; portsDiff_probeBlockA];
            portsDiff_probeBlockB_all = [portsDiff_probeBlockB_all; portsDiff_probeBlockB];
            portsDiff_noProbeBlockA_all = [portsDiff_noProbeBlockA_all; portsDiff_noProbeBlockA];
            portsDiff_noProbeBlockB_all = [portsDiff_noProbeBlockB_all; portsDiff_noProbeBlockB];
    
            idxProbe_blockA_all = [idxProbe_blockA_all; idxProbe_blockA];
            idxProbe_blockB_all = [idxProbe_blockB_all; idxProbe_blockB];
            idxNoProbe_blockA_all = [idxNoProbe_blockA_all; idxNoProbe_blockA];
            idxNoProbe_blockB_all = [idxNoProbe_blockB_all; idxNoProbe_blockB];

            cd ..
        end
        cd ..
    end

    edges = [-5.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 5.5];
    diffs = [-3, -2, -1, 0, 1, 2, 3];
    
    [portsDiff_probeBlockA_counts,~] = ...
        histcounts(portsDiff_probeBlockA_all, 'BinEdges', edges);
    [portsDiff_probeBlockB_counts,~] = ...
        histcounts(portsDiff_probeBlockB_all, 'BinEdges', edges);
    [portsDiff_noProbeBlockA_counts,~] = ...
        histcounts(portsDiff_noProbeBlockA_all, 'BinEdges', edges);
    [portsDiff_noProbeBlockB_counts,~] = ...
        histcounts(portsDiff_noProbeBlockB_all, 'BinEdges', edges);
    
    portsDiff_probeBlockA_prop = portsDiff_probeBlockA_counts / ...
        length(idxProbe_blockA_all);
    portsDiff_probeBlockB_prop = portsDiff_probeBlockB_counts / ...
        length(idxProbe_blockB_all);
    portsDiff_noProbeBlockA_prop = portsDiff_noProbeBlockA_counts / ...
        length(idxNoProbe_blockA_all);
    portsDiff_noProbeBlockB_prop = portsDiff_noProbeBlockB_counts / ...
        length(idxNoProbe_blockB_all);

    subaxis(2*length(cutoffFreq), 2, p, 'SpacingVert',0.08)
    bar(diffs, portsDiff_probeBlockA_prop, 'FaceColor', 'b')
    xlabel('Distance to target port (ports)')
    ylabel('% probe trials')
    ylim([0,1])
    title('% probe trials with target #1-3')
    
    subaxis(2*length(cutoffFreq), 2, p+1, 'SpacingVert',0.08)
    bar(diffs, portsDiff_probeBlockB_prop, 'FaceColor', 'b')
    xlabel('Distance to target port (ports)')
    ylabel('% probe trials')
    ylim([0,1])
    title('% probe trials with target #4-6')
    
    subaxis(2*length(cutoffFreq), 2, p+2, 'SpacingVert',0.08)
    bar(diffs, portsDiff_noProbeBlockA_prop, 'FaceColor', 'k')
    xlabel('Distance to target port (ports)')
    ylabel('% normal trials')
    ylim([0,1])
    title('% normal trials with target #1-3')
    
    subaxis(2*length(cutoffFreq), 2, p+3, 'SpacingVert',0.08)
    bar(diffs, portsDiff_noProbeBlockB_prop, 'FaceColor', 'k')
    xlabel('Distance to target port (ports)')
    ylabel('% normal trials')
    ylim([0,1])
    title('% normal trials with target #4-6')

    annotation('textbox', box_pos, 'String', ...
        strcat('Cutoff frequency: ', string(cutoffFreq(i)), 'kHz'), ...
        'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'EdgeColor','none');

    box_pos(2) = box_pos(2) - 0.295;
    p = p + 4;
end
        
saveas(gcf, fullfile(rootDir, 'distance_target_all'), 'fig');
saveas(gcf, fullfile(rootDir, 'distance_target_all'), 'png');
end