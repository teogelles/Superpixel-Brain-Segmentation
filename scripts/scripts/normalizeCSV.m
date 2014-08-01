% Simple function which takes the CSV files that we have generated
% from the slic features and makes them normalized for x,y, and z
% on a range from 0 to 100


function normalizeCSV(direcNum, cols)
    
    if ~exist('direcNum','var')
        direcNum = 120;
    end
    
    global debug;
    debug = false;
    
    filebase = strcat('/scratch/tgelles1/summer2014/slicExact', ...
                      num2str(direcNum), '/features/CSV/');
    savefilebase = ['/scratch/tgelles1/summer2014/slicExact', num2str(direcNum), '/features/' ...
                    'CSV_Changed/'];
    
    if ~exist(savefilebase,'dir')
        mkdir(savefilebase);
    end
    
    listing = dir(filebase);

    for i = 1:length(listing)
        file = listing(i).name;
        
        if length(file) < 4
            % not a .csv file
            ddisp('file is no good',file)
            continue
        end
        
        if ~strcmp(file(end-3:end),'.csv')
            % also not a csv
            ddisp('file is no good',file,file(end-3:end))
            continue
        end        
        
        fullfile = strcat(filebase,file);
        disp(fullfile)
        csvmat = csvread(fullfile);
        
        if ~exist('cols','var')
            cols = 1:size(csvmat,2);
        end
        csvmat = processSV(csvmat, cols);
        savefile = strcat(savefilebase,file);
        csvwrite(savefile,csvmat);
    end
end

function newSV = processSV(superVoxels,cols)
    
%   newSV = removeBackgroundSV(superVoxels);
    newSV = normalizeXYZ(superVoxels);
    % keep only wanted columns
    newSV = newSV(:,cols);
end

function newSV = normalizeXYZ(newSV)
    
    mins = [min(newSV(:,1)) min(newSV(:,2)) min(newSV(:,3))];
    maxes = [max(newSV(:,1)) max(newSV(:,2)) max(newSV(:,3))];
    for xyz = 1:3
        % normalizing x,y, and z to be from 0 to 1
        ddisp(newSV(1,xyz));
        newSV(:,xyz) = (newSV(:,xyz) - mins(xyz))/...
            (maxes(xyz) - mins(xyz));
        ddisp(newSV(1,xyz));
        % changing from 0-1 to 0-100
        % newSV(:,xyz) = 100*newSV(:,xyz);
        ddisp(newSV(1,xyz));
        % normalizing the spreads to be from 0 to 1
        newSV(:,xyz + 3) = newSV(:,xyz + 3)/(maxes(xyz) - ...
                                             mins(xyz));
        % changing spreads to also be 0 - 100
        % newSV(:,xyz + 3) = 100*newSV(:,xyz+3);
    end
end

function brainSV = removeBackgroundSV(superVoxels)
    
%threshold can be adjusted based on testing, but since we're only
%trying to remove background, it should be very low intensity
    threshold = .01;
    findByIntensity = false;
    
    if findByIntensity
        minIntensity = min(superVoxels(:,10));
        maxIntensity = max(superVoxels(:,10));
        intensityThreshold = minIntensity + (maxIntensity - minIntensity)*threshold;

        
        brainInd = find(superVoxels(:,10) > intensityThreshold);
        brainSV = superVoxels(brainInd,:);
        backgroundInd = find(superVoxels(:,10) < intensityThreshold);
        backgroundSV = superVoxels(backgroundInd,:);
        
        %creating bounding box for finding ventricals
        boxX = [min(brainSV(:,1)),max(brainSV(:,1))];
        boxY = [min(brainSV(:,2)),max(brainSV(:,2))];
        boxZ = [min(brainSV(:,3)),max(brainSV(:,3))];
        
        ventricals = find((boxX(1) < backgroundSV(:,1)) & ...
                          (boxX(2) > backgroundSV(:,1)) & ...
                          (boxY(1) < backgroundSV(:,2)) & ...
                          (boxY(2) > backgroundSV(:,2)) & ...
                          (boxZ(1) < backgroundSV(:,3)) & ...
                          (boxZ(2) > backgroundSV(:,3)));
        if size(ventricals,1) ~= 0
            numBrSV = size(brainSV,1);
            numVentSV = size(ventricals,1);
            numFeat = size(brainSV,2);
            newbrain = zeros(numBrSV + numVentSV, numFeat);
            newbrain(1:numBrSV,:) = brainSV;
            newbrain((numBrSV+1):end,:) = backgroundSV(ventricals,:);
        end
        fprintf('There are %d ventricals SVs in this brain\n', ...
                size(ventricals,1));
        fprintf('There are %d nonvent black SVs in this brain\n',...
                size(backgroundSV,1) - size(ventricals,1));
    else
        brainInd = find((superVoxels(:,7) >= threshold) | ...
                        (superVoxels(:,8) >= threshold) | ...
                        (superVoxels(:,9) >= threshold));
        brainSV = superVoxels(brainInd,:);
        % fprintf('%d of %d SV kept\n',size(brainInd,1), ...
        %         size(superVoxels,1));
    end
    
end