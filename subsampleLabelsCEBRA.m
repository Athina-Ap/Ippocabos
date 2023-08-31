function subsampleLabelsCEBRA(varargin)

p = inputParser;
addParameter(p,'root',pwd,@isstr);

parse(p,varargin{:});
root = p.Results.root;

% contents = dir(basepath);
% sessions = contents(~[contents.isdir]);
sessions = {'IZ39\Final\IZ39_220624_sess10', 'IZ39\Final\IZ39_220702_sess14', ...
    'IZ39\Final\IZ39_220714_sess18', 'IZ40\Final\IZ40_220708_sess17', ...
    'IZ40\Final\IZ40_220714_sess18', 'IZ43\Final\IZ43_220828_sess4', ...
    'IZ43\Final\IZ43_220913_sess11', 'IZ43\Final\IZ43_220919_sess14', ...
    'IZ44\Final\IZ44_220830_sess7', 'IZ44\Final\IZ44_220913_sess11', ... 
    'IZ44\Final\IZ44_220919_sess14', 'IZ47\Final\IZ47_230710_sess25'};

for s = 1:length(sessions)

    %% Load data
    cd(sessions{s})
    basepath = pwd;
    basename = basenameFromBasepath(basepath);

    % Load cebra labels
    if ~isempty(dir([basepath filesep '*.labelsCEBRA.mat']))
        file = dir([basepath filesep '*.labelsCEBRA.mat']);
        load(file.name)
    end

    % Load tracking and behavior 
    if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat']))
        file = dir([basepath filesep '*.Tracking.Behavior.mat']);
        load(file.name)
    end

    if ~isempty(dir([basepath filesep '*.TrialBehavior.Behavior.mat']))
        file = dir([basepath filesep '*.TrialBehavior.Behavior.mat']);
        load(file.name)
    end

    % Load spikes
    if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
        file = dir([basepath filesep '*.spikes.cellinfo.mat']);
        load(file.name)
    end
    
    %% Find smallest trial
    % Find indices of start of first trial and end of last trial 
    dt = mean(diff(tracking.timestamps));
    win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
    spkData = bz_SpktToSpkmat(spikes,'dt',dt,'win',win);
    
    numTrials = size(behavTrials.timestamps,1);
    
    startIdx = [];
    endIdx = [];
    for tt = 1:numTrials
        [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
        [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
    end

    % Find length of all fwd runs for tone trials
    idxFwd = find(~isnan(relDistStop));
    fwd_trial_len = [];
    startFwdIdx = [];
    for tt = 1:numTrials
        commonIdx = intersect(idxFwd, startIdx(tt):endIdx(tt));   
        if ~isempty(commonIdx)
            [fwd_trial_len(tt)] = length(commonIdx);
            [startFwdIdx(tt)] = commonIdx(1);
        else
            [fwd_trial_len(tt)] = nan;
            [startFwdIdx(tt)] = nan;
        end
    end
    
    % Find min length overall (in case it is not trial 1)
    min_trial_len = min(fwd_trial_len(fwd_trial_len~=0));

    %% Subsample labels 

    % Find new sample indices
    new_samples = [];
    for i = 1:numTrials
        if ~isnan(fwd_trial_len(i)) & fwd_trial_len(i) ~= min_trial_len
            subsample_idx = linspace(1, fwd_trial_len(i), min_trial_len);
            new_samples = [new_samples, round(interp1(startFwdIdx(i):startFwdIdx(i)+fwd_trial_len(i), subsample_idx, 'linear'))];
        elseif fwd_trial_len(i) == min_trial_len
            new_samples = [new_samples, startFwdIdx(i):startFwdIdx(i)+fwd_trial_len(i)-1];
        else
            continue
        end
    end
    
    % Get new labels
    data = spkMat(:,new_samples);
    position = y(1, new_samples);
    direction = y(2, new_samples);
    trial = trialType(new_samples);
    distance = relDistStop(new_samples);

    save([basepath filesep basename '.subLabelsCEBRA.mat'], ...
        'data', 'position', 'direction', 'trial', 'distance');

    cd(root)

end
end