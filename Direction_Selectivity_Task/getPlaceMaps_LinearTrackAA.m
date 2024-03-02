function [positions, idxL_tracked, idxR_tracked] = getTrackedPositions(behavTrials, tracking, iter, sessLen)

%% Defaults and Parms


%% Bin positions in left vs right running directions
% Left lick is followed by right run and vice versa.
% Only consider the first run away from the port, in case the animal runs
% back without triggering the opposite detector. 

dirChangeR = cell(sessLen,1);
dirChangeL = cell(sessLen,1);

% First indices of each direction 
idxR_first = [];
idxL_first = [];

% Indices of all positions considered 
idxR_tracked = [];
idxL_tracked = [];

% Position threshold: note that this is relative to the starting position,
% because sometimes the starting position is in the middle of the port
% because of filtering. 
% max ~480, min ~25
% thresL = 400; % pixels away from right port
% thresR = 105; % pixels away from left port
if iter ~= 0
    start = behavTrials.start(iter);
else
    start = string(behavTrials.start(:)');
end

for i = sessLen-1
    if strcmp(start, 'left')
        [idxR] = InIntervals(tracking.timestamps,behavTrials.timestamps(i,:));
        currentTracking = tracking.position.y(idxR);
        currentTimestamps = tracking.timestamps(idxR);
        idxR_first = [idxR_first, find(diff(idxR) ~= 0, 1, "first") + 1]; % first trial idx
        thresL = currentTracking(1) - 80; % relative threshold
        dirCond = find(diff(sign(tracking.position.vy(idxR))) ~= 0); % direction condition
        posCond = find(tracking.position.y(idxR) < thresL); % crossing position threshold condition
        dirChangeR{i} = intersect(dirCond, posCond);
        
        if ~isempty(dirChangeR{i})
            positions.right{i} = [currentTimestamps(1:dirChangeR{i}(1)) currentTracking(1:dirChangeR{i}(1))];
            idxR_tracked = [idxR_tracked, idxR_first(i):idxR_first(i)+dirChangeR{i}(1)];
        else
            positions.right{i} = [tracking.timestamps(idxR) tracking.position.y(idxR)];
            idxR_tracked = [idxR_tracked, find(idxR==1)'];
        end
        
        [idxL] = InIntervals(tracking.timestamps,[behavTrials.timestamps(i,2) behavTrials.timestamps(i+1,1)]);
        currentTracking = tracking.position.y(idxL);
        currentTimestamps = tracking.timestamps(idxL);
        idxL_first = [idxL_first, find(diff(idxL) ~= 0, 1, "first") + 1]; % first trial idx
        thresR = currentTracking(1) + 80;
        dirCond = find(diff(sign(tracking.position.vy(idxL))) ~= 0);
        posCond = find(tracking.position.y(idxL) > thresR);
        dirChangeL{i} = intersect(dirCond, posCond);

        if ~isempty(dirChangeL{i})
            positions.left{i} = [currentTimestamps(1:dirChangeL{i}(1)) currentTracking(1:dirChangeL{i}(1))];
            idxL_tracked = [idxL_tracked, idxL_first(i):idxL_first(i)+dirChangeL{i}(1)];
        else
            positions.left{i} = [tracking.timestamps(idxL) tracking.position.y(idxL)];
            idxL_tracked = [idxL_tracked, find(idxL==1)'];
        end

    else
        [idxL] = InIntervals(tracking.timestamps,behavTrials.timestamps(i,:));
        currentTracking = tracking.position.y(idxL);
        currentTimestamps = tracking.timestamps(idxL);
        idxL_first = [idxL_first, find(diff(idxL) ~= 0, 1, "first") + 1]; % first trial idx
        thresR = currentTracking(1) + 80;
        dirCond = find(diff(sign(tracking.position.vy(idxL))) ~= 0); 
        posCond = find(tracking.position.y(idxL) > thresR); 
        dirChangeL{i} = intersect(dirCond, posCond);

        if ~isempty(dirChangeL{i})
            positions.left{i} = [currentTimestamps(1:dirChangeL{i}(1)) currentTracking(1:dirChangeL{i}(1))];
            idxL_tracked = [idxL_tracked, idxL_first(i):idxL_first(i)+dirChangeL{i}(1)];
        else
            positions.left{i} = [tracking.timestamps(idxL) tracking.position.y(idxL)];
            idxL_tracked = [idxL_tracked, find(idxL==1)'];
        end

        [idxR] = InIntervals(tracking.timestamps,[behavTrials.timestamps(i,2) behavTrials.timestamps(i+1,1)]);
        currentTracking = tracking.position.y(idxR);
        currentTimestamps = tracking.timestamps(idxR);
        idxR_first = [idxR_first, find(diff(idxR) ~= 0, 1, "first") + 1]; % first trial idx
        thresL = currentTracking(1) - 80;
        dirCond = find(diff(sign(tracking.position.vy(idxR))) ~= 0); 
        posCond = find(tracking.position.y(idxR) < thresL); 
        dirChangeR{i} = intersect(dirCond, posCond);

        if ~isempty(dirChangeR{i})
            positions.right{i} = [currentTimestamps(1:dirChangeR{i}(1)) currentTracking(1:dirChangeR{i}(1))];
            idxR_tracked = [idxR_tracked, idxR_first(i):idxR_first(i)+dirChangeR{i}(1)];
        else
            positions.right{i} = [tracking.timestamps(idxR) tracking.position.y(idxR)];
            idxR_tracked = [idxR_tracked, find(idxR==1)'];
        end    
    end
end

% Fix last trial
if strcmp(start, 'left')
    % last lick is right
    [idxR] = InIntervals(tracking.timestamps,behavTrials.timestamps(end,:));
    currentTracking = tracking.position.y(idxR);
    currentTimestamps = tracking.timestamps(idxR);
    idxR_first = [idxR_first, find(diff(idxR) ~= 0, 1, "first") + 1]; % first trial idx
    thresL = currentTracking(1) - 80;
    dirCond = find(diff(sign(tracking.position.vy(idxR))) ~= 0); 
    posCond = find(tracking.position.y(idxR) < thresL); 
    dirChangeR{end} = intersect(dirCond, posCond);
    
    if ~isempty(dirChangeR{end})
        positions.right{end} = [currentTimestamps(1:dirChangeR{end}(1)) currentTracking(1:dirChangeR{end}(1))];
        idxR_tracked = [idxR_tracked, idxR_first(end):idxR_first(end)+dirChangeR{end}(1)];
    else
        positions.right{end} = [tracking.timestamps(idxR) tracking.position.y(idxR)];
        idxR_tracked = [idxR_tracked, find(idxR==1)'];
    end    
else
    % last lick is left
    [idxL] = InIntervals(tracking.timestamps,behavTrials.timestamps(end,:));
    currentTracking = tracking.position.y(idxL);
    currentTimestamps = tracking.timestamps(idxL);
    idxL_first = [idxL_first, find(diff(idxL) ~= 0, 1, "first") + 1]; % first trial idx
    thresR = currentTracking(1) + 80;
    dirCond = find(diff(sign(tracking.position.vy(idxL))) ~= 0); 
    posCond = find(tracking.position.y(idxL) > thresR); 
    dirChangeL{end} = intersect(dirCond, posCond);
    
    if ~isempty(dirChangeL{end})
        positions.left{end} = [currentTimestamps(1:dirChangeL{end}(1)) currentTracking(1:dirChangeL{end}(1))];
        idxL_tracked = [idxL_tracked, idxL_first(end):idxL_first(end)+dirChangeL{end}(1)];
    else
        positions.left{end} = [tracking.timestamps(idxL) tracking.position.y(idxL)];
        idxL_tracked = [idxL_tracked, find(idxL==1)'];
    end
end

end