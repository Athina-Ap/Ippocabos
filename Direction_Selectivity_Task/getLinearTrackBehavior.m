function [behavTrials] = getLinearTrackBehavior(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',false,@islogical)
addParameter(p,'forceRun',true,@islogical)

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;

basepath = pwd;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.TrialBehavior.Events.mat'])) && ~forceRun
    disp('Trial behavior already detected! Loading file.');
    file = dir([basepath filesep '*.TrialBehavior.Events.mat']);
    load(file.name);
    return
end

if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) 
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
end

%% Get digital inputs
if exist('settings.xml')
    delete 'settings.xml'
end
disp('Loading digital In...');
digitalIn = bz_getDigitalIn_AA;
if isempty(digitalIn)
    toneBehav = [];
    return
end

%% Read digital inputs
% Left solenoid is digital input 16, right solenoid is digital input 10 
digitalIn.timestampsOn{10} = digitalIn.timestampsOn{10}(digitalIn.dur{10}>0.01);
digitalIn.timestampsOn{16} = digitalIn.timestampsOn{16}(digitalIn.dur{16}>0.01);

% TODO: Only consider behavior timestamps after the tracking has started 
if (digitalIn.timestampsOn{16}(1) < tracking.timestamps(1)) | ...
        (digitalIn.timestampsOn{10}(1) < tracking.timestamps(1))

    disp('Oops, the mouse licked before the tracking started. It is being corrected...')
    [lickR] = find(digitalIn.timestampsOn{10} > tracking.timestamps(1), 1, 'first');
    [lickL] = find(digitalIn.timestampsOn{16} > tracking.timestamps(1), 1, 'first');
    
    if digitalIn.timestampsOn{16}(lickL) < digitalIn.timestampsOn{10}(lickR)
        behavTrials.start = 'left';
    else
        behavTrials.start = 'right';
    end

    if length(digitalIn.timestampsOn{16}(lickL:end)) == length(digitalIn.timestampsOn{10}(lickR:end))
        behavTrials.timestamps = [digitalIn.timestampsOn{16}(lickL:end) digitalIn.timestampsOn{10}(lickR:end)];   
    elseif length(digitalIn.timestampsOn{16}(lickL:end)) > length(digitalIn.timestampsOn{10}(lickR:end))
        % start was left 
        behavTrials.timestamps = [digitalIn.timestampsOn{16}(lickL:end-1) digitalIn.timestampsOn{10}(lickR:end)];
    elseif length(digitalIn.timestampsOn{16}(lickL:end)) < length(digitalIn.timestampsOn{10}(lickR:end))
        % start was right 
        behavTrials.timestamps = [digitalIn.timestampsOn{10}(lickR:end-1) digitalIn.timestampsOn{16}(lickL:end)];    
    end

else
    if digitalIn.timestampsOn{16}(1) < digitalIn.timestampsOn{10}(1)
        behavTrials.start = 'left';
    else
        behavTrials.start = 'right';
    end
    
    if length(digitalIn.timestampsOn{16}) == length(digitalIn.timestampsOn{10})
        behavTrials.timestamps = [digitalIn.timestampsOn{16}(1:end) digitalIn.timestampsOn{10}(1:end)];   
    elseif length(digitalIn.timestampsOn{16}) > length(digitalIn.timestampsOn{10})
        % start was left 
        behavTrials.timestamps = [digitalIn.timestampsOn{16}(1:end-1) digitalIn.timestampsOn{10}(1:end)];
    elseif length(digitalIn.timestampsOn{16}) < length(digitalIn.timestampsOn{10})
        % start was right 
        behavTrials.timestamps = [digitalIn.timestampsOn{10}(1:end-1) digitalIn.timestampsOn{16}(1:end)];    
    end
end

behavTrials.numTrials = size(behavTrials.timestamps,1)*2;

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehavior.Events.mat'],'behavTrials');
end

end