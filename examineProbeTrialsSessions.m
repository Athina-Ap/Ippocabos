% This function examines the behavior of the mice during the probe trials.
% At the moment, the probe trials are trials where the tone stops playing
% when it reaches a certain frequency. This function combines data across
% sessions for each animal. 

function [probeBehavior] = examineProbeTrialsSessions(varargin)

p = inputParser;
addParameter(p,'rootDir',pwd,@isstr); % Folder with all sessions
parse(p,varargin{:});
rootDir = p.Results.rootDir;

files = dir(rootDir);
sessions = files([files.isdir]);
sessions = sessions(~ismember({sessions(:).name},{'.','..'}));

cutoffFreq = cell2mat(extractBetween(rootDir,'_','kHz')); % kHz

% Pre-define variables for all sessions 
toneGain_sessions = [];

idxProbe_sessions = [];
idxNoProbe_sessions = [];
idxCorrectProbe_sessions = [];
idxCorrectNoProbe_sessions = [];
idxIncorrectProbe_sessions = [];
idxIncorrectNoProbe_sessions = [];

portsDiff_probeBlockA_sessions = [];
portsDiff_probeBlockB_sessions = [];
portsDiff_noProbeBlockA_sessions = [];
portsDiff_noProbeBlockB_sessions = [];

idxProbe_blockA_sessions = [];
idxProbe_blockB_sessions = [];
idxNoProbe_blockA_sessions = [];
idxNoProbe_blockB_sessions = [];

for s = 1:length(sessions)
    basepath = fullfile(sessions(s).folder, sessions(s).name);

    % Load the behavior file 
    if ~isempty(dir([basepath filesep '*TrialBehavior.Events.mat'])) 
        disp('Behavior already detected! Loading file.');
        file = dir([basepath filesep '*TrialBehavior.Events.mat']);
        load(fullfile(basepath, file(1).name));
    end
    
    %% Determine licking location and correct vs incorrect trials
    idxProbe = find(behavTrials.probeTrial==1);
    idxNoProbe = find(behavTrials.probeTrial==0);
    toneGain = behavTrials.toneGain + 1;
    lickLoc = behavTrials.lickLoc + 1; 
    
    idxProbe_sessions = [idxProbe_sessions; idxProbe];
    idxNoProbe_sessions = [idxNoProbe_sessions; idxNoProbe];
    toneGain_sessions = [toneGain_sessions; toneGain];

    idxCorrectProbe = find(behavTrials.correct==1 & behavTrials.toneTrial==0 ...
        & behavTrials.linTrial==0 & behavTrials.probeTrial==1);
    idxCorrectNoProbe = find(behavTrials.correct==1 & behavTrials.toneTrial==0 ...
        & behavTrials.linTrial==0 & behavTrials.probeTrial==0);
    
    idxIncorrectProbe = find(behavTrials.correct==0 & behavTrials.toneTrial==0 ...
        & behavTrials.linTrial==0 & behavTrials.probeTrial==1);
     idxIncorrectNoProbe = find(behavTrials.correct==0 & behavTrials.toneTrial==0 ...
        & behavTrials.linTrial==0 & behavTrials.probeTrial==0);
    
    idxCorrectProbe_sessions = [idxCorrectProbe_sessions; idxCorrectProbe];
    idxCorrectNoProbe_sessions = [idxCorrectNoProbe_sessions; idxCorrectNoProbe];
    idxIncorrectProbe_sessions = [idxIncorrectProbe_sessions; idxIncorrectProbe];
    idxIncorrectNoProbe_sessions = [idxIncorrectNoProbe_sessions; idxIncorrectNoProbe];
    
    %% Calculate the distance to target based on the location of cutoffFreq
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

    % Update variables for all sessions
    portsDiff_probeBlockA_sessions = [portsDiff_probeBlockA_sessions; ...
        portsDiff_probeBlockA];
    portsDiff_probeBlockB_sessions = [portsDiff_probeBlockB_sessions; ...
        portsDiff_probeBlockB];
    portsDiff_noProbeBlockA_sessions = [portsDiff_noProbeBlockA_sessions; ...
        portsDiff_noProbeBlockA];
    portsDiff_noProbeBlockB_sessions = [portsDiff_noProbeBlockB_sessions; ...
        portsDiff_noProbeBlockB];
    
    idxProbe_blockA_sessions = [idxProbe_blockA_sessions; idxProbe_blockA];
    idxProbe_blockB_sessions = [idxProbe_blockB_sessions; idxProbe_blockB];
    idxNoProbe_blockA_sessions = [idxNoProbe_blockA_sessions; idxNoProbe_blockA];
    idxNoProbe_blockB_sessions = [idxNoProbe_blockB_sessions; idxNoProbe_blockB];
end

%% Visualize distribution of licks (cor/incor) and targets for probe trials
edges = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5];
[licksProbe_sessions,~] = histcounts(toneGain_sessions(idxProbe_sessions), ...
    'BinEdges', edges);
[licksNoProbe_sessions,~] = histcounts(toneGain_sessions(idxNoProbe_sessions), ...
    'BinEdges', edges);
[licksCorrectProbe_sessions,~] = histcounts(toneGain_sessions ...
    (idxCorrectProbe_sessions), 'BinEdges', edges);
[licksCorrectNoProbe_sessions,~] = histcounts(toneGain_sessions ...
    (idxCorrectNoProbe_sessions), 'BinEdges', edges);
[licksIncorrectProbe_sessions,~] = histcounts(toneGain_sessions ...
    (idxIncorrectProbe_sessions), 'BinEdges', edges);
[licksIncorrectNoProbe_sessions,~] = histcounts(toneGain_sessions ...
    (idxIncorrectNoProbe_sessions), 'BinEdges', edges);

figure('WindowState','maximized');

subplot(2,3,1)
bar(1:6, licksProbe_sessions);
xlabel('Port')
ylabel('Number of probe trials')
title('Number of probe trials for each target port')

subplot(2,3,2)
bar(1:6, licksCorrectProbe_sessions);
xlabel('Port')
ylabel('Number of correct probe trials')
title('Number of correct probe trials for each target port')

subplot(2,3,3)
bar(1:6, licksIncorrectProbe_sessions);
xlabel('Port')
ylabel('Number of incorrect probe trials')
title('Number of incorrect probe trials for each target port')

subplot(2,3,4)
bar(1:6, licksNoProbe_sessions);
xlabel('Port')
ylabel('Number of non-probe trials')
title('Number of non-probe trials for each target port')

subplot(2,3,5)
bar(1:6, licksCorrectNoProbe_sessions);
xlabel('Port')
ylabel('Number of correct non-probe trials')
title('Number of correct non-probe trials for each target port')

subplot(2,3,6)
bar(1:6, licksIncorrectNoProbe_sessions);
xlabel('Port')
ylabel('Number of incorrect non-probe trials')
title('Number of incorrect non-probe trials for each target port')

saveas(gcf, fullfile(rootDir, strcat('probe_corr_incorr_sessions_',cutoffFreq,'kHz')), 'fig');
saveas(gcf, fullfile(rootDir, strcat('probe_corr_incorr_sessions_',cutoffFreq,'kHz')), 'png');

%% Plot the distribution of distance to target for all sessions
edges = [-5.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 5.5];
diffs = [-3, -2, -1, 0, 1, 2, 3];

[portsDiff_probeBlockA_counts_sessions,~] = ...
    histcounts(portsDiff_probeBlockA_sessions, 'BinEdges', edges);
[portsDiff_probeBlockB_counts_sessions,~] = ...
    histcounts(portsDiff_probeBlockB_sessions, 'BinEdges', edges);
[portsDiff_noProbeBlockA_counts_sessions,~] = ...
    histcounts(portsDiff_noProbeBlockA_sessions, 'BinEdges', edges);
[portsDiff_noProbeBlockB_counts_sessions,~] = ...
    histcounts(portsDiff_noProbeBlockB_sessions, 'BinEdges', edges);

portsDiff_probeBlockA_prop_sessions = portsDiff_probeBlockA_counts_sessions / ...
    length(idxProbe_blockA_sessions);
portsDiff_probeBlockB_prop_sessions = portsDiff_probeBlockB_counts_sessions / ...
    length(idxProbe_blockB_sessions);
portsDiff_noProbeBlockA_prop_sessions = portsDiff_noProbeBlockA_counts_sessions / ...
    length(idxNoProbe_blockA_sessions);
portsDiff_noProbeBlockB_prop_sessions = portsDiff_noProbeBlockB_counts_sessions / ...
    length(idxNoProbe_blockB_sessions);

figure('WindowState','maximized');
subplot(2,2,1)
bar(diffs, portsDiff_probeBlockA_prop_sessions)
xlabel('Distance to target port (ports)')
ylabel('% probe trials (targets #1-3)')
title(['Sessions: Proportion of probe trials with target #1-3 for each distance ' ...
    'to target port'])

subplot(2,2,2)
bar(diffs, portsDiff_probeBlockB_prop_sessions)
xlabel('Distance to target port (ports)')
ylabel('% probe trials (targets #4-6)')
title(['Sessions: Proportion of probe trials with target #4-6 for each distance ' ...
    'to target port'])

subplot(2,2,3)
bar(diffs, portsDiff_noProbeBlockA_prop_sessions, 'FaceColor', 'r')
xlabel('Distance to target port (ports)')
ylabel('% non-probe trials (targets #1-3)')
title(['Sessions: Proportion of non-probe trials with target #1-3 for each distance ' ...
    'to target port'])

subplot(2,2,4)
bar(diffs, portsDiff_noProbeBlockB_prop_sessions, 'FaceColor', 'r')
xlabel('Distance to target port (ports)')
ylabel('% non-probe trials (targets #4-6)')
title(['Sessions: Proportion of non-probe trials with target #4-6 for each distance ' ...
    'to target port'])

saveas(gcf, fullfile(rootDir, strcat('distance_target_sessions_',cutoffFreq,'kHz')), 'fig');
saveas(gcf, fullfile(rootDir, strcat('distance_target_sessions_',cutoffFreq,'kHz')), 'png');

end