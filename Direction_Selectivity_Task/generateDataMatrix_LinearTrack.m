function generateDataMatrix_LinearTrack(varargin)

basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

%% Load data 
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Tracking already detected! Loading file.');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(fullfile(basepath, file(1).name));
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(fullfile(basepath, file(1).name));
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(fullfile(basepath, file.name));
end
%% (1) Generate spike matrix 
% To do that, we use the dt as defined by the tracking.timestamps. 

dt = mean(diff(tracking.timestamps));
win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
spkData = bz_SpktToSpkmat(spikes,'dt',dt,'win',win);

spkMat = spkData.data';
timestamps = spkData.timestamps';

%% (2) Generate variables 
idx_behav = InIntervals(tracking.timestamps, ...
    [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);

y = tracking.position.y(idx_behav); 
% Correct mismatches in the sizes
if length(y) ~= size(spkMat(1,:))
    idx_behav(find(idx_behav==1, 1, 'first') - 1) = 1;
end
y = tracking.position.y(idx_behav);

% left = 0, right = 1
direction = ones(size(y)); 
direction(tracking.position.vy(idx_behav) < 0) = 0; 

start = behavTrials.start;
event.times{1} = behavTrials.timestamps(:,1);
event.times{2} = behavTrials.timestamps(:,2);
event.UID(1) = 1;
event.UID(2) = 2;
eventMat = bz_SpktToSpkmat(event,'dt',dt,'win',win);

if strcmp(start, 'left')
    licksL = eventMat.data(:,1);
    licksR = eventMat.data(:,2);
else
    licksL = eventMat.data(:,2);
    licksR = eventMat.data(:,1);
end

%% (3) Save data for CEBRA
save(fullfile(basepath, strcat(basename, '.labelsCEBRA.mat')), ...
        'y', 'direction', 'licksR', 'licksL', 'spkMat', 'start');

end