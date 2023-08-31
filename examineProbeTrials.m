% This function examines the behavior of the mice during the probe trials.
% At the moment, the probe trials are trials where the tone stops playing
% when it reaches a certain frequency. 

function [probeBehavior] = examineProbeTrials(varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
parse(p,varargin{:});
basepath = p.Results.basepath;

% Load the behavior file 
if ~isempty(dir([basepath filesep '*TrialBehavior.Events.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Events.mat']);
    load(fullfile(basepath, file(1).name));
end

cutoffFreq = cell2mat(extractBetween(basepath,'_','kHz')); % kHz - check spreadsheet for this

% Check the session had probe trials
if ~isfield(behavTrials, 'probeTrial')
    disp('This session had no probe trials.')
    return
end
%% Visualize distribution of licks (cor/incor) and targets for probe trials
idxProbe = find(behavTrials.probeTrial==1);
idxNoProbe = find(behavTrials.probeTrial==0);
toneGain = behavTrials.toneGain + 1;
lickLoc = behavTrials.lickLoc + 1; 

figure('WindowState','maximized');

edges = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5];
% Distribution of sampling ports for probe trials 
[licksProbe,~] = histcounts(toneGain(idxProbe), 'BinEdges', edges);
[licksNoProbe,~] = histcounts(toneGain(idxNoProbe), 'BinEdges', edges);

subplot(2,3,1)
bar(1:6, licksProbe);
xlabel('Port')
ylabel('Number of probe trials')
title('Number of probe trials for each target port')

subplot(2,3,4)
bar(1:6, licksNoProbe);
xlabel('Port')
ylabel('Number of no-probe trials')
title('Number of no-probe trials for each target port')

% Distribution of correct probe trials 
idxCorrectProbe = find(behavTrials.correct==1 & behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==1);
idxCorrectNoProbe = find(behavTrials.correct==1 & behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==0);

[licksCorrectProbe,~] = histcounts(toneGain(idxCorrectProbe), 'BinEdges', edges);
[licksCorrectNoProbe,~] = histcounts(toneGain(idxCorrectNoProbe), 'BinEdges', edges);

subplot(2,3,2)
bar(1:6, licksCorrectProbe);
xlabel('Port')
ylabel('Number of correct probe trials')
title('Number of correct probe trials for each target port')

subplot(2,3,5)
bar(1:6, licksCorrectNoProbe);
xlabel('Port')
ylabel('Number of correct no-probe trials')
title('Number of correct no-probe trials for each target port')


% Distribution of incorrect probe trials 
idxIncorrectProbe = find(behavTrials.correct==0 & behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==1);
idxIncorrectNoProbe = find(behavTrials.correct==0 & behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==0);

[licksIncorrectProbe,~] = histcounts(toneGain(idxIncorrectProbe), 'BinEdges', edges);
[licksIncorrectNoProbe,~] = histcounts(toneGain(idxIncorrectNoProbe), 'BinEdges', edges);

subplot(2,3,3)
bar(1:6, licksIncorrectProbe);
xlabel('Port')
ylabel('Number of incorrect probe trials')
title('Number of incorrect probe trials for each target port')

subplot(2,3,6)
bar(1:6, licksIncorrectNoProbe);
xlabel('Port')
ylabel('Number of incorrect no-probe trials')
title('Number of incorrect no-probe trials for each target port')

saveas(gcf, fullfile(basepath, strcat('licks_targets_',string(cutoffFreq),'kHz')), 'fig');
saveas(gcf, fullfile(basepath, strcat('licks_targets_',string(cutoffFreq),'kHz')), 'png');

%% Calculate accuracy for probe and no-probe trials
accuracy_noProbe = length(idxCorrectNoProbe) / sum([length(idxCorrectNoProbe), ...
    length(idxIncorrectNoProbe)]);
accuracy_Probe = length(idxCorrectProbe) / sum([length(idxCorrectProbe), ...
    length(idxIncorrectProbe)]);

if accuracy_noProbe >= 0.65
    disp(['Accuracy for non-probe trials is greater than 65%. Probe trials ' ...
        'can be considered.'])
else
    disp(['Accuracy for non-probe trials is less than 65%. Please do not ' ...
        'consider the probe trials in these sessions.'])
end

probeTrials = {};
probeTrials.accuracy_noProbe = accuracy_noProbe;
probeTrials.accuracy_Probe = accuracy_Probe;

save(fullfile(basepath, strcat(extractAfter(basepath, 'kHz\'), ...
    '.probeTrials.mat')), 'probeTrials'); 

%% Plot the distribution of distance to lick based on the location of cutoffFreq
% Distance = target - lickloc so negative distance means the lick is after
% the target.

% Distance to target for probe trials. Block A = 1-3, Block B = 4-6
idxProbe_blockA = find(behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==1 & ...
    (toneGain == 1 | toneGain == 2 | toneGain == 3));

idxProbe_blockB = find(behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==1 & ...
    (toneGain == 4 | toneGain == 5 | toneGain == 6));

portsDiff_probeBlockA = zeros(length(idxProbe_blockA), 1);
portsDiff_probeBlockB = zeros(length(idxProbe_blockB), 1);

for i = 1:length(idxProbe_blockA)
    portsDiff_probeBlockA(i) = toneGain(idxProbe_blockA(i)) - ...
        lickLoc(idxProbe_blockA(i));
end

for i = 1:length(idxProbe_blockB)
    portsDiff_probeBlockB(i) = toneGain(idxProbe_blockB(i)) - ...
        lickLoc(idxProbe_blockB(i));
end

edges = [-5.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 5.5];
[portsDiff_probeBlockA_counts,~] = histcounts(portsDiff_probeBlockA, 'BinEdges', edges);
[portsDiff_probeBlockB_counts,~] = histcounts(portsDiff_probeBlockB, 'BinEdges', edges);

portsDiff_probeBlockA_prop = portsDiff_probeBlockA_counts / length(idxProbe_blockA);
portsDiff_probeBlockB_prop = portsDiff_probeBlockB_counts / length(idxProbe_blockB);

% diffs = ceil(edges);
% diffs(end) = [];
diffs = [-3, -2, -1, 0, 1, 2, 3];

figure('WindowState','maximized');
subplot(2,2,1)
bar(diffs, portsDiff_probeBlockA_prop)
xlabel('Distance to target port (ports)')
ylabel('% probe trials (targets #1-3)')
title(['Proportion of probe trials with target #1-3 for each distance ' ...
    'to target port'])

subplot(2,2,2)
bar(diffs, portsDiff_probeBlockB_prop)
xlabel('Distance to target port (ports)')
ylabel('% probe trials (targets #4-6)')
title(['Proportion of probe trials with target #4-6 for each distance ' ...
    'to target port'])

% Distance to target for non-probe trials. Block A = 1-3, Block B = 4-6
idxNoProbe_blockA = find(behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==0 & ...
    (toneGain == 1 | toneGain == 2 | toneGain == 3));

idxNoProbe_blockB = find(behavTrials.toneTrial==0 ...
    & behavTrials.linTrial==0 & behavTrials.probeTrial==0 & ...
    (toneGain == 4 | toneGain == 5 | toneGain == 6));

portsDiff_noProbeBlockA = zeros(length(idxNoProbe_blockA), 1);
portsDiff_noProbeBlockB = zeros(length(idxNoProbe_blockB), 1);

for i = 1:length(idxNoProbe_blockA)
    portsDiff_noProbeBlockA(i) = toneGain(idxNoProbe_blockA(i)) - ...
        lickLoc(idxNoProbe_blockA(i));
end

for i = 1:length(idxNoProbe_blockB)
    portsDiff_noProbeBlockB(i) = toneGain(idxNoProbe_blockB(i)) - ...
        lickLoc(idxNoProbe_blockB(i));
end

[portsDiff_noProbeBlockA_counts,~] = histcounts(portsDiff_noProbeBlockA, 'BinEdges', edges);
[portsDiff_noProbeBlockB_counts,~] = histcounts(portsDiff_noProbeBlockB, 'BinEdges', edges);

portsDiff_noProbeBlockA_prop = portsDiff_noProbeBlockA_counts / length(idxNoProbe_blockA);
portsDiff_noProbeBlockB_prop = portsDiff_noProbeBlockB_counts / length(idxNoProbe_blockB);

subplot(2,2,3)
bar(diffs, portsDiff_noProbeBlockA_prop, 'FaceColor', 'r')
xlabel('Distance to target port (ports)')
ylabel('% non-probe trials (targets #1-3)')
title(['Proportion of non-probe trials with target #1-3 for each distance ' ...
    'to target port'])

subplot(2,2,4)
bar(diffs, portsDiff_noProbeBlockB_prop, 'FaceColor', 'r')
xlabel('Distance to target port (ports)')
ylabel('% non-probe trials (targets #4-6)')
title(['Proportion of non-probe trials with target #4-6 for each distance ' ...
    'to target port'])

saveas(gcf, fullfile(basepath, strcat('distance_target_',string(cutoffFreq),'kHz')), 'fig');
saveas(gcf, fullfile(basepath, strcat('distance_target_',string(cutoffFreq),'kHz')), 'png');

end