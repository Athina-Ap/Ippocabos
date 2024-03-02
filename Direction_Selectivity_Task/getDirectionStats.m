%% (1) getDirectionStats (place cells + stability among directions)
function [directionStats] = getDirectionStats(labels, numcells, placeFieldStats)

% Get cells with place field(s)
for d = 1:length(labels)
    placeCells{d} = [];
    for i = 1:numcells 
        if ~isnan(placeFieldStats.mapStats{i}{d}.x)
            placeCells{d} = [placeCells{d}, i];
        end
    end
    directionStats.placeCells.(labels{d}) = placeCells{d};
    directionStats.props.(labels{d}) = length(placeCells{d}) / numcells;
end

for d = 1:length(labels)
    for i = 1:length(placeCells{d})
        directionStats.mapStats{d}{i} = placeFieldStats.mapStats{placeCells{d}(i)}{d};
        directionStats.mapStats{d}{i}.UID = placeCells{d}(i);
    end
end

% Get cells with fields in two directions
cells2fields = intersect(placeCells{1}, placeCells{2});
directionStats.cells2fields.UID = cells2fields;

% Find overlapping fields
for unit = 1:length(cells2fields)
    numFields = [];
    for d = 1:length(labels)
        numFields = [numFields, size(placeFieldStats.mapStats{cells2fields(unit)}{d}.field, 2)];
    end

    directionStats.cells2fields.sameField{unit} = zeros(1,max(numFields));
    for d = 1
        od = 3-d; % opposite direction

        for ii = 1:numFields(d)
            for jj = 1:numFields(od)
                if length(intersect(find(placeFieldStats.mapStats{cells2fields(unit)}{od}.field(:,jj)), ...
                        find(placeFieldStats.mapStats{cells2fields(unit)}{d}.field(:,ii)))) >= ...
                        (0.8 * min(length(find(placeFieldStats.mapStats{cells2fields(unit)}{od}.field(:,jj))), ...
                        length(find(placeFieldStats.mapStats{cells2fields(unit)}{d}.field(:,ii)))))
                    
                    directionStats.cells2fields.sameField{unit}(ii) = 1;
                end
            end
        end
    end
end

directionStats.cells2fields.samePerc = length(find(cell2mat(directionStats.cells2fields.sameField) == 1)) / ...
    length(cell2mat(directionStats.cells2fields.sameField));

end