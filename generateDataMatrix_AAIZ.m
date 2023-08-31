%% Variables of interest: 
% x position on the track
% y position on the track (circular) - tone and no-tone trials
% Absolute distance to the stop location
% Licking type: tone/choice, tone/spontaneous, no-tone/choice,
%   no-tone/spontaneous, home 
% Licking events 
% Theta LFP 
%
% Spike counts 30 ms
% 
% Spike matrix

function [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix_AAIZ(varargin)

%% Defaults and Parms - Change for each session
p = inputParser;
addParameter(p,'basepath','I:\Videos\IZ43_220828_sess4',@isstr);
%addParameter(p,'basepath',pwd,@isstr);
%addParameter(p,'dt',1/30,@isnumeric); % Change this to 1/30 s

parse(p,varargin{:});
basepath = p.Results.basepath;
%dt = p.Results.dt;

%% (1) Load all the .mat files
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
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

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

% Load digitalIn
filename = [basepath filesep sessionInfo.session.name '.DigitalIn.events.mat'];
load(filename);

%% (2) Change what defines the beginning and the end of the trials. 
% Until now,it was defined based on DigitalIn {10} and {2} respectively. 
% Now, we want to define them based on the licking, as there is a ~33 ms 
% delay between the two signals. 

indStart = [];
indStart(1) = 1; 

% Note that the threshold might need to be different in some cases. Check
% that the sizes of all variables match in the end. 
for i = 2:length(digitalIn.timestampsOn{3})
    if abs(digitalIn.timestampsOn{3}(i)-digitalIn.timestampsOn{3}(i-1))>=5
        indStart = [indStart, i];
    end
end

% Remove the last element, as this would mark the end of the last trial and
% not the beginning of a new. 
indStart(end) = [];

behavTrials.timestampsNew(:,1) = digitalIn.timestampsOn{3}(indStart);

% Find the median delay between the lick and solenoid
delay = median(digitalIn.timestampsOn{10}(1:30) - ...
    behavTrials.timestampsNew(1:30));

% Adjust and shift the timestamps by the median delay between a lick and
% the solenoid
behavTrials.timestamps = behavTrials.timestamps - abs(delay);

%% (3) Generate spike matrix 
% To do that, we use the dt as defined by the tracking.timestamps. 

dtime = mean(diff(tracking.timestamps));

win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'win',win);

spkMat = spkData.data';
timestamps = spkData.timestamps';

%% (4) Generate matrix for continuous variables x, vel, velY, cyclicY, cyclicYLin.
% Note, that we later change the length of these variables to match the
% length of the behavior variables in case frames have been dropped. 

% Find indices of start of first trial and end of last trial 
[~,idxStart] = min(abs(tracking.timestamps-timestamps(1)));
[~,idxEnd] = min(abs(tracking.timestamps-timestamps(end)));

constVar.x = tracking.position.x(idxStart:idxEnd)'; % [cm]
constVar.y = tracking.position.y(idxStart:idxEnd)'; % [cm]
constVar.velY = tracking.position.vy(idxStart:idxEnd)'; % [cm/s]
constVar.hd = tracking.position.angle(idxStart:idxEnd)'; % [degrees]

% Convert the y position to cover both directions in the range [0,240].
constVar.ylong = constVar.y;
indices = find(constVar.velY < 0 & (constVar.hd > 85 | constVar.hd < -85));
constVar.ylong(indices) = ...
    240 - constVar.ylong(indices);
constVar.ylong = smooth(constVar.ylong, 20, 'rlowess');

constVar.cyclicY = mod(constVar.ylong, 240); % Tone trials
constVar.cyclicYLin = constVar.cyclicY; % No-tone trials

%% Find index in matrix to start and end points of each trial
numTrials = size(behavTrials.timestamps,1);

for tt = 1:numTrials
    [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
    [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
end

%% (5) Generate logical variables - 
% trial number, forward/reverse, toneOn/Off, trial gain, current choice, 
% correct/incorrect, past choice, past choice correct/incorrect

logVar.forward(1:length(spkData.timestamps)) = 0;
logVar.toneOn(1:length(spkData.timestamps)) = 0;
logVar.trialType(1:length(spkData.timestamps)) = 0;
logVar.currChoice(1:length(spkData.timestamps)) = 0;
logVar.currCorrect(1:length(spkData.timestamps)) = 0;
logVar.trialNum(1:length(spkData.timestamps)) = 0;

%last lin track trial 
lastLin = behavTrials.linTrial;
toneGain = (behavTrials.toneGain);
toneGain(~lastLin) = toneGain(~lastLin)+1;
toneGain(lastLin==1) = 0;
choiceCorrect = behavTrials.correct;
choiceCorrect(lastLin==1) = NaN;

for tt = 1:(length(startIdx)-1)        
    if tt<length(startIdx)-1
        logVar.toneOn(startIdx(tt):startIdx(tt+1)) = ~behavTrials.linTrial(tt);
        logVar.trialNum(startIdx(tt):startIdx(tt+1)) = tt;
        logVar.trialType(startIdx(tt):startIdx(tt+1)) = toneGain(tt);
        logVar.currChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt)+1;
        logVar.currCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt);
    elseif tt == length(startIdx)-1
        logVar.toneOn(startIdx(tt):endIdx(tt)) = ~behavTrials.linTrial(tt);
        logVar.trialNum(startIdx(tt):endIdx(tt)) = tt;
        logVar.trialType(startIdx(tt):endIdx(tt)) = toneGain(tt);
        logVar.currChoice(startIdx(tt):endIdx(tt)) = behavTrials.lickLoc(tt)+1;
        logVar.currCorrect(startIdx(tt):endIdx(tt)) = choiceCorrect(tt);
    end
end

%% (7) Now can generate the remaining continuous variables:

% a. Distance to stop location (absolute distance). It is only defined in 
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

%% (8) Generate matrix for events - licks, trialStart, trialEnd

lickport = [1 2 3 4 5 6 7];
for ll = 1:7
    event.times{ll} = digitalIn.timestampsOn{lickport(ll)+2}; 
    event.UID(ll) = ll;
end

event.times{8} = behavTrials.timestamps(:,1);
event.times{9} = behavTrials.timestamps(:,2);

choiceIncorrect = behavTrials.correct == 0 & behavTrials.linTrial == 0;
event.times{10} = behavTrials.timestamps(choiceIncorrect,2);

choiceCorrect = behavTrials.correct == 1;
event.times{11} = behavTrials.timestamps(choiceIncorrect,2);
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

% Create a vector with single licks when they first appear - logical. 
eventVar.lickEvents = zeros(size(eventVar.licksAll));
eventVar.lickEvents(singlesIdx) = 1;

%% (9) Separate the types of licks: 
% 1. tone choice lick, 2. tone spontaneous lick, 3. tone end lick, 
% 4. no-tone spontaneous lick, 5. home lick 

% 1. All tone trial choice licks
% ind1 = find(eventVar.licksPortsAllSingle == (logVar.currChoice+1) & ...
%     eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==1); 
ind1 = find(eventVar.licksPortsAllSingle == (logVar.currChoice+1) & ...
    eventVar.licksPortsAllSingle ~= 0 & ...
    eventVar.licksPortsAllSingle ~= 1 & logVar.toneOn==1); 

% currChoice goes from 1-6, while licks are from port 1 (home) to 7. add 1

% 2. All tone trial spontaneous licks
ind2 = find(eventVar.licksPortsAllSingle ~= (logVar.currChoice+1) & ...
    eventVar.licksPortsAllSingle ~= 0 & ...
    eventVar.licksPortsAllSingle ~= 1 & logVar.toneOn==1);

% 3. All no tone trial end licks
ind3 = find(eventVar.licksPortsAllSingle == 7 & logVar.toneOn==0);

% 4. All no tone trial middle port licks
ind4 = find(eventVar.licksPortsAllSingle ~= 7 & ...
    eventVar.licksPortsAllSingle ~= 1 & ...
    eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==0);

% 5. All home ports
ind5 = find(eventVar.licksPortsAllSingle == 1);

% Create a single vector with a different type of lick indicated.
eventVar.licksType = zeros(size(eventVar.licksAll));
eventVar.licksType(ind1) = 1;
eventVar.licksType(ind2) = 2;
eventVar.licksType(ind3) = 3;
eventVar.licksType(ind4) = 4;
eventVar.licksType(ind5) = 5;

% Create 5 different binary vectors for each type of lick. 
eventVar.('licksType1') = zeros(size(eventVar.licksAll));
eventVar.('licksType2') = zeros(size(eventVar.licksAll));
eventVar.('licksType3') = zeros(size(eventVar.licksAll));
eventVar.('licksType4') = zeros(size(eventVar.licksAll));
eventVar.('licksType5') = zeros(size(eventVar.licksAll));

eventVar.('licksType1')(ind1) = 1;
eventVar.('licksType2')(ind2) = 1;
eventVar.('licksType3')(ind3) = 1;
eventVar.('licksType4')(ind4) = 1;
eventVar.('licksType5')(ind5) = 1;

%% (10) Generate LFP data 
cd 'IZ44_220830_sess7'

lfp = bz_GetLFP(69,'noprompts',true,'interval',win);
signal.data = bz_Filter(double(lfp.data),'filter','butter','passband',[6 12],'order',3);
signal.timestamps = lfp.timestamps;
signal.samplingRate = lfp.samplingRate;

downFactor = 20;
lfpdown1 = bz_DownsampleLFP(signal, downFactor);

[~,matchIdx] = min(abs(lfpdown1.timestamps-timestamps));

lfpdown = lfpdown1.data(matchIdx);
constVar.lfp = lfpdown';

cd ..
%% Plot distributions and time series of the variables and save data. 
% Variables of interest are: eventVar.lickType, eventVar.licks, constVar.cyclicYLin, 
% constVar.cyclicY, constVar.distStop, constVar.x, lfp (theta)

all_variables.x = constVar.x;
all_variables.cyclicYLin = constVar.cyclicYLin;
all_variables.cyclicY = constVar.cyclicY;
all_variables.distStop = constVar.distStop;
all_variables.licks = eventVar.lickEvents;
all_variables.('licksType1') = eventVar.('licksType1');
all_variables.('licksType2') = eventVar.('licksType2');
all_variables.('licksType3') = eventVar.('licksType3');
all_variables.('licksType4') = eventVar.('licksType4');
all_variables.('licksType5') = eventVar.('licksType5');
all_variables.lfp = constVar.lfp;

save(fullfile(basepath, 'sessionData.mat'),'eventVar','constVar', ...
    'timestamps','spkMat','tracking','all_variables'); 

% Create plots
% fields = fieldnames(all_variables);
% 
% figure;
% 
% nc = 3;
% nr = ceil(length(fields)/3);
% 
% for i = 1:length(fields)
%     subplot(nr, nc, i);
% 
%     fieldname = fields{i};
%     field = all_variables.(fieldname);
% 
%     if strcmp(fieldname, 'licksType1') || strcmp(fieldname, 'licksType2')...
%             || strcmp(fieldname, 'licksType3') || strcmp(fieldname, 'licksType4')...
%             || strcmp(fieldname, 'licksType5') || strcmp(fieldname, 'licks')
%         field = nonzeros(field);
%     end
%     
%     histogram(field);    
%     xlabel(fieldname);
%     ylabel('density');
% end
% 
% saveas(gcf, fullfile(basepath, 'Variables_distribution'), 'fig');
% saveas(gcf, fullfile(basepath, 'Variables_distribution'), 'png');
% 
% figure;
% set(gcf,'Color','w')
% nr = length(fields);
% 
% for i = 1:length(fields)
%     subplot(nr, 1, i);
% 
%     fieldname = fields{i};
%     field = all_variables.(fieldname);
% 
%     if strcmp(fieldname, 'lfp')
%         plot(timestamps(1:length(constVariables.lfp)),field)
%     else
%         plot(timestamps,field);
%     end
% 
%     if strcmp(fieldname, 'licksType1')
%         title(strcat(fieldname,': Tone choice')) 
%     elseif strcmp(fieldname, 'licksType2')
%         title(strcat(fieldname,': Tone spont.')) 
%     elseif strcmp(fieldname, 'licksType3')
%         title(strcat(fieldname,': No tone choice')) 
%     elseif strcmp(fieldname, 'licksType4')
%         title(strcat(fieldname,': No-tone spont.'))
%     elseif strcmp(fieldname, 'licksType5')
%         title(strcat(fieldname,': home')) 
%     else
%         title(fieldname)
%     end
%     ylabel(field);
%     xlim([timestamps(1) timestamps(end)])
% 
%     if i < length(fields)
%         set(gca,'xtick',[])
%     end
% 
% %     axis off
% end
% xlabel('time (s)');
% 
% saveas(gcf, fullfile(basepath, 'Variables_time_series'), 'fig');
% saveas(gcf, fullfile(basepath, 'Variables_time_series'), 'png');

%% Output for the PGAM
eventVariables.trialStart = eventVar.trialStart;
eventVariables.trialEnd = eventVar.trialEnd;    
eventVariables.('licksType1') = eventVar.('licksType1');
eventVariables.('licksType2') = eventVar.('licksType2');
eventVariables.('licksType3') = eventVar.('licksType3');
eventVariables.('licksType4') = eventVar.('licksType4');
eventVariables.('licksType5') = eventVar.('licksType5');
eventVariables.licks = eventVar.lickEvents;
constVariables.x = constVar.x;
constVariables.cyclicY = constVar.cyclicY';
constVariables.cyclicYLin = constVar.cyclicYLin';
constVariables.distStop = constVar.distStop;
constVariables.lfp = constVar.lfp;

save(fullfile(basepath, strcat(extractAfter(basepath, 'Videos\'), '.sessionDataPGAM.mat')), ...
    'eventVariables', 'constVariables','spkMat'); 

%% Output for CEBRA
x = constVar.x;

% Forward direction - 1, backward direction - 0.
direction = ones(size(constVar.y));
direction(constVar.ylong > 120) = 0;

y = zeros(2, length(constVar.y));
y(1,:) = constVar.y;
y(2,:) = direction;
distStop = constVar.distStop;
lfp = constVar.lfp;
trialStart = eventVar.trialStart;
trialEnd = eventVar.trialEnd;
licksType1 = eventVar.('licksType1');
licksType2 = eventVar.('licksType2');
licksType3 = eventVar.('licksType3');
licksType4 = eventVar.('licksType4');
licksType5 = eventVar.('licksType5');
licks = eventVar.lickEvents;

save(fullfile(basepath, strcat(extractAfter(basepath, 'Videos\'), '.labelsCEBRA.mat')), ...
    'x', 'y', 'direction', 'distStop', 'lfp', 'trialStart', 'trialEnd', ...
    'licksType1', 'licksType2', 'licksType3', 'licksType4', 'licksType5', ...
    'licks', 'idxStart', 'idxEnd', 'spkMat');

end