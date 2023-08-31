
function [tracking] = Pos2Tracking_LinearTrack(aviFile,varargin)
% Get position tracking
%
% USAGE
%
%   [behavior] = Pos2Tracking(varagin)
%
% INPUTS
%   aviFile     Avi format video. If not provided, look for it in the
%   basepath f
%   It requires an avi format video and a digitalin.dat file in the
%   basepath folder.
% 
%   (OPTIONAL)
%   basePath       -(default: pwd) basePath for the recording file, in
%                   buzcode format.
%   roiTracking    - 2 x R, where 1C is x and 2C is y. By default it
%                   considers the whole video. With the option 'manual' allows to draw
%                   a ROI.
%   roiLED         - 2 x R, where 1C is x and 2C is y.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveFrames     - Creates mat file containin all frames (default false).
%   forceReload    - default false.
%   bazlerTTL      - Rx1 with timestamps from bazler ttl pulses. If not
%                   provided try to extract ttl pulses from digitalIn
%                   channel 1. If it fails, gives video time.
%   saveMat        - default true
%   artifactThreshold - max allow movements per frame (in cm, default 3).
%                   Disabled if not convFact is provided.
% 
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   x               - x position in cm/ normalize
%   y               - y position in cm/ normalize
%   timestamps      - in seconds, if Bazler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%   samplingRate    - in Hz
%   averageFrame    - 
%   
%
%   Athina Apostolelli 2023
%   The code is adapted from LED2Tracking.m (Manuel Valero) and
%   Pos2Tracking.m (Ipshita Zutshi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'fs',30,@isnumeric);
addParameter(p,'artifactThreshold',10,@isnumeric);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiTracking','manual',@ischar);
addParameter(p,'roiLED','manual',@ischar);
addParameter(p,'ledSync',[],@ismatrix);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveFrames',true,@islogical)
addParameter(p,'verbose',false,@islogical);
addParameter(p,'thresh',0.1,@isnumeric)
addParameter(p,'bazlerTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
% addParameter(p,'RGBChannel',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
artifactThreshold = p.Results.artifactThreshold;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
convFact = p.Results.convFact;
forceReload = p.Results.forceReload;
ledSync = p.Results.ledSync;
saveFrames = p.Results.saveFrames;
verbose = p.Results.verbose;
thresh = p.Results.thresh;
bazlerTtl = p.Results.bazlerTTL;
saveMat = p.Results.saveMat;
% RGBChannel = p.Results.RGBChannel;

%% Deal with inputs
% if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
%     disp('Trajectory already detected! Loading file.');
%     file = dir([basepath filesep '*Tracking.Behavior.mat']);
%     load(file.name);
%     return
% end

if ~exist('aviFile') || isempty(aviFile)
    if ~isempty(dir([basepath filesep 'test*avi']))
        aviFile = dir([basepath filesep 'test*avi']); 
        aviFile = erase(aviFile.name,'.avi');
    else
        warning('No video file!!');
        tracking = [];
        return
    end
end

% Load the red channel from the video
if ~exist([basepath filesep aviFile '.mat'],'file') || forceReload
    disp('Get average frame...');
    videoObj = VideoReader([aviFile '.avi']);   % get video
    numFrames = get(videoObj, 'NumFrames');
    frames.r = [];
    tic
    f = waitbar(0,'Getting frames...');
    for ii = 1:numFrames-1
        waitbar(ii/numFrames,f)
        temp_frames = read(videoObj,ii);        % get all frames
        frames.r(:,:,ii) = temp_frames(:,:,1);
%         frames.r = cat(3,frames.r,squeeze(temp_frames(:,:,1,:)));          % collect red, add new nested variable for blue and green if needed
    end
    temp_frames = read(videoObj,Inf);
    frames.r(:,:,end) = temp_frames(:,:,1);
    close(f)
    toc
    
    if saveFrames
        disp('Saving frames...');
        save([basepath filesep aviFile '.mat'],'frames','-v7.3');
    end
else
    disp('Loading frames from mat file...');
    load([basepath filesep aviFile '.mat'],'frames');
end

% get average frame
average_frame = mean(frames.r,3);                                          % get average frames

%% Figure out LED for sync
% if ledSync
%     
%     if ~exist([basepath filesep aviFile '_red.mat'],'file') || forceReload
%         disp('Get average frame...');
%         videoObj = VideoReader([aviFile '.avi']);   % get video
%         numFrames = get(videoObj, 'NumFrames');
%         frames = [];
%         tic
%         f = waitbar(0,'Getting frames...');
%         for ii = 1:numFrames
%             waitbar(ii/numFrames,f)
%             temp_frames = read(videoObj,ii);        % get all frames
%             frames_red(:,:,ii) = temp_frames(:,:,1);      % convert to grayscale
%         end
%         close(f)
%         toc
%         
%         if saveFrames
%             disp('Saving frames...');
%             save([basepath filesep aviFile '_red.mat'],'frames_red','-v7.3');
%         end
%     else
%         disp('Loading frames from mat file...');
%         load([basepath filesep aviFile '.mat'],'frames_red');
%     end
% end

% deal with the ROI for the tracking 
cd(basepath); cd ..; upBasepath = pwd; pwd; cd ..; up_upBasepath = pwd;cd(basepath);
if exist([basepath filesep 'roiTracking.mat'],'file')
    load([basepath filesep 'roiTracking.mat'],'roiTracking');
elseif exist([upBasepath filesep 'roiTracking.mat'],'file')
    load([upBasepath filesep 'roiTracking.mat'],'roiTracking');
    disp('ROI tracking from sessions folder... copying locally...');
    save([basepath filesep 'roiTracking.mat'],'roiTracking');
elseif exist([up_upBasepath filesep 'roiTracking.mat'],'file')
    load([up_upBasepath filesep 'roiTracking.mat'],'roiTracking');
    disp('ROI tracking from master folder... copying locally...');
    save([basepath filesep 'roiTracking.mat'],'roiTracking');
elseif isempty(roiTracking)
    roiTracking = [1 1 size(frames.r,2) size(frames.r,2) 1 ;1 size(frames.r,1) size(frames.r,1) 1 1 ]';
elseif ischar(roiTracking) && strcmpi(roiTracking,'manual')
    disp('Draw ROI for tracking...');
    h1 = figure;
    imagesc(average_frame); colormap hot; caxis([0 4*mean(average_frame(:))]);
    title('Draw ROI for tracking... (include all posible locations)','FontWeight','normal');
    roi = drawpolygon;
    roiTracking = [roi.Position; roi.Position(1,:)];
    save([basepath filesep 'roiTracking.mat'],'roiTracking');
    close(h1);
end

% deal with the ROI for the LED
if exist([basepath filesep 'roiLED.mat'],'file')
    load([basepath filesep 'roiLED.mat'],'roiLED');
elseif exist([upBasepath filesep 'roiLED.mat'],'file')
    load([upBasepath filesep 'roiLED.mat'],'roiLED');
    disp('ROI LED from master folder... copying locally...');
    save([basepath filesep 'roiLED.mat'],'roiLED');
elseif ischar(roiLED) && strcmpi(roiLED,'manual')
    disp('Draw ROI for LED...');
    h1 = figure;
%     imshow(average_frame);
    imagesc(average_frame); colormap hot;
    roi = drawpolygon;
    roiLED = [roi.Position; roi.Position(1,:)];
    save([basepath filesep 'roiLED.mat'],'roiLED');
    close(h1);
end

if isempty(convFact)                             % if convFact not provided, normalize to 1 along the longest axis
    convFact = 1/max([size(frames,1) size(frames,2)]);
    artifactThreshold = Inf;
end
xMaze = [0 size(frames,2) * convFact];
yMaze = [0 size(frames,1) * convFact];


%% DETECT LED POSITION
% bw = uint8(poly2mask(roiTracking(:,1),roiTracking(:,2),size(frames.r,1),size(frames.r,2)));
bw = poly2mask(roiTracking(:,1),roiTracking(:,2),size(frames.r,1),size(frames.r,2));
disp('Detect LED position...');
thr_fr = thresh * 255;
if ~verbose
    tic
    f = waitbar(0,'Detecting LED position...');
    for ii = 1:size(frames.r,3)
        waitbar(ii/size(frames.r,3),f)
        fr = frames.r(:,:,ii).*bw;
        bin_fr = imbinarize(double(fr),thr_fr); %
        bin_fr = bwareafilt(bin_fr,[10 300]);
        stats_fr = regionprops(bin_fr);
        maxBlob = find([stats_fr.Area]== max([stats_fr.Area]),1);
        if ~isempty(maxBlob)
            sz_fr(ii) = stats_fr(maxBlob).Area;
            Rr_x(ii) = stats_fr(maxBlob).Centroid(1);
            Rr_y(ii) = stats_fr(maxBlob).Centroid(2);
        else
            sz_fr(ii) = NaN;
            Rr_x(ii) = NaN;
            Rr_y(ii) = NaN;
        end
    end
    close(f)
    toc
else
    h1 = figure;
    hold on
    tic
    for ii = 1:size(frames.r,3)
        fr = frames.r(:,:,ii).*bw;
        bin_fr = imbinarize(double(fr),thr_fr); %
        bin_fr = bwareafilt(bin_fr,[10 300]);
        stats_fr = regionprops(bin_fr);
        maxBlob = find([stats_fr.Area]== max([stats_fr.Area]),1);
        if ~isempty(maxBlob)
            sz_fr(ii) = stats_fr(maxBlob).Area;
            Rr_x(ii) = stats_fr(maxBlob).Centroid(1);
            Rr_y(ii) = stats_fr(maxBlob).Centroid(2);
            cla
            imagesc(fr)
            plot(Rr_x(ii),Rr_y(ii),'or')
            drawnow;
        else
            sz_fr(ii) = NaN;
            Rr_x(ii) = NaN;
            Rr_y(ii) = NaN;
        end
    end
    toc
    close(h1);
end

pos = [Rr_x; Rr_y]';


%% postprocessing of LED position 
pos = pos * convFact;                                   % cm or normalized
art = find(sum(abs(diff(pos))>artifactThreshold,2))+1;  % remove artefacs as movement > 10cm/frame
pos(art,:) = NaN;

xt = linspace(0,size(pos,1)/fs,size(pos,1));            % kalman filter
[t,x,y,vx,vy,ax,ay] = trajectory_kalman_filter(pos(:,1)',pos(:,2)',xt,0);
art = find(sum(abs(diff([x y]))>artifactThreshold,2))+1;
art = [art - 2 art - 1 art art + 1 art + 2];
x(art(:)) = NaN; y(art(:)) = NaN;
F = fillmissing([x y],'linear');
x = F(:,1); y = F(:,2);

xDim = [min(roiTracking(:,1)), max(roiTracking(:,1))];
yDim = [min(roiTracking(:,2)), max(roiTracking(:,2))];
h2 = figure;
hold on
scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
caxis([t(1) t(end)])
xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
xlim(xDim); ylim(yDim);
mkdir('Behavior');
saveas(h2,'Behavior\trajectory.png');
if ~verbose
    close(h2);
end

%% detect LED pulses for sync
% if ~isempty(roiLED)
%     disp('Detect LED for sync...');
%     bwLED = uint8(poly2mask(roiLED(:,1),roiLED(:,2),size(frames,1),size(frames,2)));
%     parfor ii = 1:size(frames,3)
%         fr = double(frames(:,:,ii).*bwLED);
%         fr(fr==0) = NaN;
%         sync(ii) = nanmedian(fr(:)); % Manu uses nansum
%     end
% 
%     sync = sync.^2;
%     syncBin = (sync>median(sync)); % binarize signal
%     locsA = find(diff(syncBin)==1)/fs; % start of pulses
%     locsB = find(diff(syncBin)==-1)/fs; % end of pulses
%     pul = locsA(1:min([length(locsA) length(locsB)]));
%     for ii = 1 : size(pul,2) % pair begining and end of the pulse
%         if sum(locsB > pul(1,ii)) > 0
%             pul(2,ii) =  locsB(find(locsB - pul(1,ii) ==...
%                 min(locsB(locsB > pul(1,ii)) - pul(1,ii))));
%         else
%             pul(2,ii) = nan;
%         end
%     end
%     % if a jump happened only for 1-2 frames, ignore
%     a = pul(2,:)-pul(1,:);
% else
%     sync = []; pul = [];
% end

%% Get basler TTL
% digital_in legend: 1. Basler, 2. maze LED, 3. Left Alternation, 4.Right
% Alternation, 5. Home Delay, 6. Is alternation forzed?
cd ..; upPath = pwd; cd(basepath);
if isempty(bazlerTtl)
%     if ~isempty(dir([upPath filesep '*DigitalIn.events.mat']))
%     digitalIn = bz_getDigitalIn;
    if ~isempty(dir('*digitalIn.events.mat'))
        disp('DigitalIn already detected. Loading file.')
%         file = dir([upPath filesep '*DigitalIn.events.mat']);
        file = dir('*digitalIn.events.mat');
        load(file.name);
%         cd(basepath);
        bazlerTtl = digitalIn.timestampsOn{1};
    else
        disp('No DigitalIn file. Trying to make it now...')
        digitalIn = bz_getDigitalIn_AA;
        file = dir('*digitalIn.events.mat');
        load(file.name);
        bazlerTtl = digitalIn.timestampsOn{1};
    end
end

tracking.framesDropped = length(bazlerTtl) - length(x);

% match basler frames with ttl pulses
if length(bazlerTtl) == length(x)
    disp('Number of frames match!!');
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) > length(x) + 5 && ~isempty(dir([basepath filesep 'test*csv'])) %% if a csv file exists
    
    txtFile = dir([basepath filesep 'test*csv']);
    fid = fopen(txtFile.name, 'r');
    behavInfo = textscan(fid,'%s','delimiter','\n');
    TS = cellfun(@(x)strsplit(x,'_'),behavInfo{1},'UniformOutput',false);       
    TS = vertcat(TS{:});
    TS = str2double(TS);
    TS = round(TS-TS(1));
    if length(TS) ~= length(x)
        disp('Yikes!! CSV file does not match');
    end
    %Determine how the dropped frames align
    tsBas = round((bazlerTtl-bazlerTtl(1))*1000);    
    Xcorr(1:length(tsBas)) = nan;
    Ycorr(1:length(tsBas)) = nan;
    
    for tt = 1:length(TS)
        [~, idx] = min(abs(tsBas-TS(tt)));
        Xcorr(idx) = x(tt);
        Ycorr(idx) = y(tt);
    end
    
    fprintf('%3.i frames were dropped, fixed using csv file \n',...
        length(bazlerTtl) - length(x));
    x = Xcorr;
    y = Ycorr;
    
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) <= length(x) + 30 * 1 
    fprintf('%3.i frames were dropped, probably at the end of the recording. Skipping... \n',...
        length(bazlerTtl) - length(x));
    bazlerTtl = bazlerTtl(1:length(x));
elseif length(bazlerTtl) > length(x) && length(bazlerTtl) > length(x) + 30 * 1 
    fprintf('%3.i frames were dropped!! Examine the recording carefully!... \n',...
        length(bazlerTtl) - length(x));
elseif length(bazlerTtl) < length(x) 
    fprintf('%3.i video frames without TTL... Was the recording switched off before the camera?. Skipping... \n',...
        length(x) - length(bazlerTtl));
    x = x(1:length(bazlerTtl));
    y = y(1:length(bazlerTtl));
    vx = vx(1:length(bazlerTtl));
    vy = vy(1:length(bazlerTtl));
    ax = ax(1:length(bazlerTtl));
    ay = ay(1:length(bazlerTtl)); 
elseif isempty(bazlerTtl)
    bazlerTtl = xt;
end

[~,fbasename,~]=fileparts(pwd);

% Get instantaneous speed
y = smooth(y, 20, 'rlowess');
[~,~,~,vx,vy,~,~] = KalmanVel(x,y,bazlerTtl,2);
v = sqrt(vx.^2+vy.^2);

tracking.position.x =x;
tracking.position.y = y;
tracking.position.vx =vx;
tracking.position.vy = vy;
tracking.position.v =v;

tracking.description = '';
tracking.timestamps = bazlerTtl;
tracking.originalTimestamps = [];
tracking.folder = fbasename;
tracking.samplingRate = fs;
tracking.avFrame.r = average_frame;
tracking.avFrame.xSize = xMaze;
tracking.avFrame.ySize = yMaze;
tracking.roi.roiTracking = roiTracking;
tracking.roi.roiLED = roiLED;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

% maze is 48 x 67 (exterior wall to exterior wall). Virtual maze is 580 x
% 420 pixels. Conv factor is ~ 0.1143 - 0.1155 cm/pix 
