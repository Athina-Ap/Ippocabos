function [tracking] = getLinearTrackTracking(varargin)

% Gets position tracking for each sub-session and concatenate all of them so they are 
% aligned with LFP and spikes. Default is using the analog input from Bonsai detected position. 
%
%% USAGE
%
%   [tracking] = getToneTracking(varargin)
%
% INPUTS
%   basePath        -(default: pwd) basePath for the recording file, in buzcode format:
%   analogInputPos  - Analog channel that has tracking information.
%   analogInputTone - Analog input that has tone information.
%   fs              - sampling rate for behavior. default 1250
%   trackImgLength  - Distance of track in the video file, default, 410.
%   trackLength     - Actual length of the track (in cm)
%   freqRange       - Frequency range of the Tone, default 1000 to 22000
%   saveMat         - default true
%   forceReload     - default false
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   position.x               - x position in cm/ normalize
%   position.tone
%   timestamps      - in seconds, if Basler ttl detected, sync by them

%   HISTORY:
%     - Ipshita Zutshi 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'analogInputPos',2,@isnumeric);
addParameter(p,'analogInputTone',3,@isnumeric);
addParameter(p,'fs',150,@isnumeric);
addParameter(p,'trackLength',112,@isnumeric);
addParameter(p,'trackImgLength',410,@isnumeric);
addParameter(p,'freqRange',[1000 22000],@isnumeric);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',true,@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
analogInputPos = p.Results.analogInputPos;
analogInputTone = p.Results.analogInputTone;
fs = p.Results.fs;
trackLength = p.Results.trackLength;
trackImgLength = p.Results.trackImgLength;
freqRange = p.Results.freqRange;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;

%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

%% Get tracking and behavior (no subfolders)
cd(basepath);
% 
% cd ..; upPath = pwd; cd(basepath);
% sessPattern = extractBetween(basepath, strcat(upPath,'\'), 'sess');
% pattern = [basepath, filesep, sessPattern, '*'];
% contents = dir(strcat(pattern{:}));
% for i = 1:length(contents)
%     if contents(i).isdir == 1
%         cd(contents(i).name)
%         tracking = Pos2Tracking_LinearTrack; % computing trajectory
%         behavTrials = getLinearTrackBehavior('forceRun',true);
%     end
% end
% 
% cd(basepath)

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'test*']))   
            cd([basepath filesep MergePoints.foldernames{ii}]); 
            fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
            tempTracking{count}= Pos2Tracking_LinearTrack; % computing trajectory
            tempBehav{count} = getLinearTrackBehavior('forceRun',true);
            trackFolder(count) = ii; 
            count = count + 1;
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quiting...');
end

%% Concatenate and sync timestamps
ts = []; subSessions = []; maskSessions = [];
tsBehav= []; maskSessionsBehav = [];
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    for ii = 1:length(trackFolder)
        if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
            sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
            maskSessions = [maskSessions; ones(size(sumTs))*ii];
            ts = [ts; sumTs];
            
            sumTs = tempBehav{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            maskSessionsBehav = [maskSessionsBehav; ones(size(sumTs))*ii];
            tsBehav = [tsBehav; sumTs];
        else
            error('Folders name does not match!!');
        end
    end
else
    warning('No MergePoints file found. Concatenating timestamps...');
    for ii = 1:length(trackFolder)
        sumTs = max(ts)+ tempTracking{ii}.timestamps;
        subSessions = [subSessions; [sumTs(1) sumTs(end)]];
        ts = [ts; sumTs];
    end
end

% Concatenating tracking fields...
x = []; y = []; vx = []; vy = []; v = []; folder = []; samplingRate = []; description = [];framesDropped = [];
for ii = 1:size(tempTracking,2) 
    x = [x; tempTracking{ii}.position.x];     
    y = [y; tempTracking{ii}.position.y]; 
    
    vx = [vx; tempTracking{ii}.position.vx];     
    vy = [vy; tempTracking{ii}.position.vy]; 
    v = [v; tempTracking{ii}.position.v];     
    
    folder{ii} = tempTracking{ii}.folder; 
    samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
    description{ii} = tempTracking{ii}.description;  
    framesDropped{ii} = tempTracking{ii}.framesDropped;  
end

tracking.position.x = x;
tracking.position.y = y;
tracking.position.vx = vx;
tracking.position.vy = vy;
tracking.position.v = v;

tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = ts;
tracking.framesDropped = framesDropped;
tracking.events.subSessions =  subSessions;
tracking.events.subSessionsMask = maskSessions;

% behavTrials = tempBehav{1};
behavTrials.start = [string(tempBehav{1}.start), string(tempBehav{2}.start)];
behavTrials.numTrials = [tempBehav{1}.numTrials, tempBehav{2}.numTrials];
behavTrials.timestamps = tsBehav;
% behavTrials.start = temp

%% save tracking 
[sessionInfo] = bz_getSessionInfo(basepath,'noPrompts', true,'saveMat',true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
    save([basepath filesep sessionInfo.FileName '.TrialBehavior.Behavior.mat'],'behavTrials');
end

end

