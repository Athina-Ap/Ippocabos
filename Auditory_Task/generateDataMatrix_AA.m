%% Variables of interest: 
% x position on the track
% y position on the track (circular) - tone and no-tone trials
% Absolute distance to the stop location
% Licking type: tone/choice, tone/spontaneous, no-tone/spontaneous, home 
% Theta LFP 
%
% Spike counts 33.3 ms
% 
% Spike matrix

function [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix_AA(varargin)

addpath(genpath('F:\github\buzcode'));

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'dt',1/30,@isnumeric); % Change this to 1/30 s

parse(p,varargin{:});
basepath = p.Results.basepath;
dt = p.Results.dt;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

% Load digitalIn
filename = [basepath filesep sessionInfo.session.name '.DigitalIn.events.mat'];
load(filename);

%% Change what defines the beginning and the end of the trials. 
% Until now,it was defined based on DigitalIn {10} and {2} respectively. 
% Now, we want to define them based on the licking, as there is a ~33 ms 
% delay between the two signals. 

indStart = [];
indStart(1) = 1; 

for i = 2:length(digitalIn.timestampsOn{3})
    if abs(digitalIn.timestampsOn{3}(i)-digitalIn.timestampsOn{3}(i-1))>=10
        indStart = [indStart, i];
    end
end

% Remove the last element, as this would mark the end of the last trial and
% not the beginning of a new. 
indStart(end) = [];

behavTrials.timestampsNew(:,1) = digitalIn.timestampsOn{3}(indStart);

delay = median(digitalIn.timestampsOn{10}(1:30) - ...
    behavTrials.timestampsNew(1:30));

behavTrials.timestamps = behavTrials.timestamps - delay;

%% First generate spike matrix 
% To do that, we use the dt as defined by the tracking.timestamps. 

dtime = mean(diff(tracking.timestamps));

win = [behavTrials.timestamps(1,1) behavTrials.timestamps((end-1),2)];
spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'win',win);

spkMat = spkData.data';
timestamps = spkData.timestamps';

%% Next, generate matrix for constant variables x, y, vel, velY, cyclicY.
% Note, that we later change the length of these variables to match the
% length of the behavior variables in case frames have been dropped. 

% Find indices of start of first trial and end of last trial 
[~,idxStart] = min(abs(tracking.timestamps-win(1)));
[~,idxEnd] = min(abs(tracking.timestamps-win(2)));

constVar.x = tracking.position.x(idxStart:idxEnd)'; % [cm]
constVar.y = tracking.position.y(idxStart:idxEnd)'; % [cm]
constVar.velY = tracking.position.vy(idxStart:idxEnd)';

% Convert the y position to cover both directions in the range [0,240].
constVar.ylong = constVar.y;
constVar.ylong(constVar.velY < -1) = 240 - constVar.ylong(constVar.velY < -1);
constVar.cyclicY = mod(constVar.ylong, 240); 
constVar.cyclicYLin = constVar.cyclicY;

%% Find index in matrix to start and end points of each trial
numTrials = size(behavTrials.timestamps,1);

for tt = 1:numTrials
    [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
    [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
end

%% Find indices of missing frames and set the frames to NaN.
% startTimestamps = timestamps(startIdx);
% [~,lickIdx] = min(abs(startTimestamps-tracking.timestamps(idxStart:idxEnd))); 
% tracking.timestampsNew = tracking.timestamps';
% 
% % Used to shift all indices by the number of new frames added.
% newNum = 0;
% 
% for t = 1:length(startIdx)-1
%     trackNum = lickIdx(t+1)-lickIdx(t);
%     behavNum = startIdx(t+1)-startIdx(t); 
% 
%     if trackNum ~= behavNum 
%         disp(['There are ', num2str(abs(trackNum-behavNum)), ...
%             ' frames missing from trial ', num2str(t)]);
% 
%         new = nan(numel(abs(trackNum-behavNum)));
% 
%         tracking.timestampsNew = [tracking.timestampsNew(...
%             1:idxStart+lickIdx(t+1)+newNum-1), new(1,:), ...
%             tracking.timestampsNew(idxStart+lickIdx(t+1)+newNum:end)];    
% 
%         newNum = newNum + numel(abs(trackNum-behavNum));
%     end
% end
% 
% trackNumEnd = lickIdx(end)-lickIdx(end-1);
% behavNumEnd = startIdx(end)-startIdx(end-1);
% 
% if trackNumEnd ~= behavNumEnd
%     disp('There are frames missing from the last trial.');
% 
%     new = nan(numel(abs(trackNumEnd-behavNumEnd)));
% 
%     tracking.timestampsNew = [tracking.timestampsNew(...
%         1:idxStart+lickIdx(end)+newNum-1), new(1,:), ...
%         tracking.timestampsNew(idxStart+lickIdx(end)+newNum:end)];    
% 
%     newNum = newNum + numel(abs(trackNumEnd-behavNumEnd));
% end
% 
% tracking.timestampsNew = tracking.timestampsNew';

% %% Change the const variables using the new timestamps. 
% % Find indices of start of first trial and end of last trial 
% idxStartNew = startIdx(1);
% idxEndNew = endIdx(end);
% 
% constVar.y = constVar.yOld;
% constVar.x = constVar.xOld;
% constVar.velY = constVar.velYOld;
% 
% empty = find(isnan(tracking.timestampsNew(idxStartNew:idxEndNew)));
% added = 0;

% Fill the tracking arrays with the median values of the frames around them
% so that the tracking matches the behavior. 
% for i = 1:length(empty)
%     constVar.y = [constVar.y(1:empty(i)+added-1), ...
%         median([constVar.yOld(empty(i)-1),constVar.yOld(empty(i))]), ...
%         constVar.y(empty(i)+added:end)];
% 
%     constVar.x = [constVar.x(1:empty(i)+added-1), ...
%         median([constVar.xOld(empty(i)-1),constVar.xOld(empty(i))]), ...
%         constVar.x(empty(i)+added:end)];
% 
%     constVar.velY = [constVar.velY(1:empty(i)+added-1), ...
%         median([constVar.velYOld(empty(i)-1),constVar.velYOld(empty(i))]), ...
%         constVar.velY(empty(i)+added:end)];
% 
%     added = added + 1;
% end 

%% Generate logical variables - 
% trial number, forward/reverse, toneOn/Off, trial gain, current choice, 
% correct/incorrect, past choice, past choice correct/incorrect

logVar.forward(1:length(spkData.timestamps)) = 0;
logVar.toneOn(1:length(spkData.timestamps)) = 0;

%last lin track trial 
lastLin = behavTrials.linTrial;
toneGain = (behavTrials.toneGain);
toneGain(~lastLin) = toneGain(~lastLin)+1;
toneGain(lastLin==1) = 0;
choiceCorrect = behavTrials.correct;
choiceCorrect(lastLin==1) = NaN;

for tt = 1:(length(startIdx)-1)
    logVar.trialNum(startIdx(tt):startIdx(tt+1)) = tt;
    logVar.toneOn(startIdx(tt):endIdx(tt)) = ~behavTrials.linTrial(tt);
    logVar.trialType(startIdx(tt):startIdx(tt+1)) = toneGain(tt);
    logVar.currChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt)+1;
    logVar.currCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt);
    if tt==1
        logVar.prevChoice(startIdx(tt):startIdx(tt+1)) = NaN;
        logVar.prevCorrect(startIdx(tt):startIdx(tt+1)) = NaN; 
    else
        logVar.prevChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt-1)+1;
        logVar.prevCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt-1);    
    end
end

%% Now can generate the remaining continuous variables:
% 1. Frequency - transformed to Hz

gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];

constVar.freq(1:length(spkData.timestamps)) = nan;

freqExp = log10(22000/1000);
for ii = 1:length(logVar.trialType)
    if logVar.trialType(ii)>0 && logVar.toneOn(ii)>0
        freq = (constVar.y(ii)*gain(logVar.trialType(ii)))/120;
         constVar.freq(ii) = ((1000*(10.^(freqExp*freq)))/22000)*120;
    end
end

% 2. Distance to stop location (absolute distance). It is only defined in 
% the boundaries of the trials and only for tone trials.  
port_locations = [11.6, 32.27, 55.53, 79.62, 102.79, 120];
constVar.distStop = NaN(size(logVar.currChoice));

for i = 1:numTrials
    if logVar.trialType(startIdx(i):endIdx(i)) > 0
        constVar.distStop(startIdx(i):endIdx(i)) =  ...
        port_locations(logVar.currChoice(startIdx(i):endIdx(i))) - ...
        constVar.y(startIdx(i):endIdx(i));
    end
end

correct = find(constVar.distStop < 0 & ~isnan(constVar.distStop));
constVar.distStop(correct) = 0;

constVar.cyclicYLin(logVar.trialType>0) = nan; % No-tone trials
constVar.cyclicY(logVar.trialType==0) = nan; % Tone trials

%% Generate matrix for events - licks, trialStart, trialEnd

lickport = [1 2 3 4 5 6 7];
for ll = 1:7
    event.times{ll} = digitalIn.timestampsOn{lickport(ll)+2}; 
    event.UID(ll) = ll;
end

event.times{8} = behavTrials.timestamps((1:(end-1)),1);
event.times{9} = behavTrials.timestamps((1:(end-1)),2);

choiceIncorrect = behavTrials.correct == 0 & behavTrials.linTrial == 0;
event.times{10} = behavTrials.timestamps(choiceIncorrect(1:(end-1)),2);

choiceCorrect = behavTrials.correct == 1;
event.times{11} = behavTrials.timestamps(choiceCorrect(1:(end-1)),2);
event.UID(8:11)  = [8 9 10 11];

eventMat = bz_SpktToSpkmat(event,'dt',dtime,'win',win);

% Add the first and last trials
eventMat.data(1,8) = 1;
eventMat.data(end,9) = 1;

% Store other trial relevant information.
eventVar.trialStart = eventMat.data(:,8)';
eventVar.trialEnd = eventMat.data(:,9)';
eventVar.incorrect = eventMat.data(:,10)';
eventVar.correct = eventMat.data(:,11)';

% Note that these licks now are only the ones within the boundaries of the
% first and last trial (win) - logical.
eventVar.licks = eventMat.data(:,1:7)';

% Create a single vector with all licking events - logical.
eventVar.licksAll = eventMat.data(:,1)' | eventMat.data(:,2)' | ...
    eventMat.data(:,3)' | eventMat.data(:,4)' | eventMat.data(:,5)' | ...
    eventMat.data(:,6)' | eventMat.data(:,7)';

% Create a vector with licking events according to port.
eventVar.licksPorts = eventMat.data(:,1:7)';
for i = 1:7
    indices = find(eventVar.licks(i,:) == 1);
    eventVar.licksPorts(i,indices) = i;
end

eventVar.licksPortsAll = zeros(size(eventVar.licksAll));
for i = 1:7
    indices = find(eventVar.licks(i,:) == 1);
    eventVar.licksPortsAll(indices) = ...
        eventVar.licksPorts(i,indices);
end

% Create a vector with the id of the lick port only at the index where it
% first appears anew. 
eventVar.licksPortsAllSingle = zeros(size(eventVar.licksAll));

licksIdx = find(eventVar.licksPortsAll ~= 0);
difference = diff(eventVar.licksPortsAll(licksIdx));
singlesIdx = licksIdx(find(difference ~=0) + 1);

eventVar.licksPortsAllSingle(singlesIdx) = eventVar.licksPortsAll(singlesIdx);
eventVar.licksPortsAllSingle(1) = 1;

% Separate the types of licks: 1. tone choice lick, 2. tone spontaneous 
% lick, 3. no-tone spontaneous lick, 4. home lick 

idxNoToneTrials = find(toneGain == 0);
idxToneTrials = find(toneGain ~= 0);

ind1 = [];
ind2 = [];
ind3 = [];

start_indices = find(eventVar.trialStart == 1);
end_indices = find(eventVar.trialEnd == 1);

tone_end = end_indices(idxToneTrials(end));
noTone_start = start_indices(idxNoToneTrials(find(diff(idxNoToneTrials)>1)+1));

toneTrials = any(idxToneTrials(:) == find(start_indices)); 

for j = 1:(length(idxToneTrials) + length(idxNoToneTrials))
    % Tone trials
    if toneTrials(j) == 1
        if j == idxToneTrials(end)
            break
        end
        ind1 = [ind1, end_indices(j)];

        indices = end_indices(j) + find(eventVar.licksPortsAllSingle(...
            (end_indices(j)+1):(start_indices(j+1)-1)) > 1);
        ind2 = [ind2, indices];
    
    % No-tone trials
    else
        if j == idxNoToneTrials(end) | j == idxToneTrials(1)
            break
        end
        indices = end_indices(j) + find(eventVar.licksPortsAllSingle(...
            (end_indices(j)+1):(start_indices(j+1)-1)) > 1);
        ind3 = [ind3, indices];
    end
end

ind2 = [ind2, tone_end + find(eventVar.licksPortsAllSingle(...
    tone_end+1:noTone_start-1) > 1)];

ind3 = [ind3, end_indices(end) + find(eventVar.licksPortsAllSingle(...
    end_indices(end)+1:end) > 1)];

ind4 = find(eventVar.licksPortsAllSingle(start_indices(1):end_indices(end)) == 1);

eventVar.licksType = zeros(size(eventVar.licksAll));
eventVar.licksType(ind1) = 1;
eventVar.licksType(ind2) = 2;
eventVar.licksType(ind3) = 3;
eventVar.licksType(ind4) = 4;

%% Generate LFP data 
lfp = bz_GetLFP(69,'basepath','I:\Videos\IZ43_220828_sess4','noprompts',true);
signal = bz_Filter(double(lfp.data),'filter','butter','passband',[6 12],'order',3);
downFactor = floor(length(signal)/length(timestamps));
lfpdown = downsample(signal,downFactor);
constVar.lfp = lfpdown(1:length(timestamps))';

%% Plot distributions and time series of the variables. 
% Variables of interest are: eventVar.lickType, constVar.cyclicYLin, 
% constVar.cyclicY, constVar.distStop, constVar.x, lfp (theta)
all_variables = struct('x', [], 'cyclicYLin', [], 'cyclicY', [], 'distStop', ...
    [], 'licksType', []);
all_variables.x = constVar.x;
all_variables.cyclicYLin = constVar.cyclicYLin;
all_variables.cyclicY = constVar.cyclicY;
all_variables.distStop = constVar.distStop;
all_variables.licksType = eventVar.licksType;
all_variables.lfp = constVar.lfp;

save(fullfile(basepath, 'sessionData.mat'),'eventVar','constVar','timestamps','spkMat',...
    'tracking','all_variables'); 

fields = fieldnames(all_variables);

figure;

nc = 3;
nr = ceil(length(fields)/3);
np = 1;

for i = 1:length(fields)
    subplot(nr, nc, np);

    fieldname = fields{i};
    field = all_variables.(fieldname);

    if strcmp(fieldname, 'licksType')
        field = nonzeros(field);
    end
    
    histogram(field);    
    xlabel(fieldname);
    ylabel('density');

    np = np + 1;
end

saveas(gcf, fullfile(basepath, 'Variables_distribution'), 'fig');
saveas(gcf, fullfile(basepath, 'Variables_distribution'), 'png');

figure;
nc = 3;
nr = ceil(length(fields)/3);
np = 1;

for i = 1:length(fields)
    subplot(nr, nc, np);

    fieldname = fields{i};
    field = all_variables.(fieldname);

    if strcmp(fieldname, 'licksType')
        bar(field);
    else
        plot(field);
    end

    xlabel('time (ms)');
    ylabel(field);

    np = np + 1;
end

saveas(gcf, fullfile(basepath, 'Variables_time_series'), 'fig');
saveas(gcf, fullfile(basepath, 'Variables_time_series'), 'png');

end