%% getSLICFeatures
% Authors: Teo Gelles & Andrew Gilchrist-Scott
%
% This file contains the code for obtaining a list of statistical
% features about a brain given various sources of information about
% the brain.  This program also writes the features to a file

function avgEntropy = getSLICFeatures(im,labels,tissues, ...
                                          centerInfo,cropOffset,filename,id)
    % getSLICFeatures - Get the list of features for the brain
    %
    % @param im - Matrix representation of an MRI of the brain
    % @param labels - Matrix representation of the SLIC labels of
    % the brain
    % @param tissues - Matrix representation of the 3 main tissue
    % types of the brain (Grey Matter, White Matter, and Cerebral
    % Spinal Fluid)
    % @param centerInfo - Collection of data about the centers used
    % in the SLIC supervoxelation of the image.  Main items of
    % interest are the total intensity of a superVoxel of given
    % center and the total number of voxels in each superVoxel
    % @param cropOffset - The dimensions by which the original
    % image was cut in order to remove the black margins around the
    % brain in cropBlack()
    % @param filename - The filename to write the features to
    % @param id - An identification number for the brain whose
    % features are being obtained
    %
    % @returns featureList - Matrix of the names of all obtained
    % feature in one column and values of the features in second
    % column.
    
    shouldPrint = true;
    fprintf('Getting SLIC Features for Patient %d\n', id);

    [avgEntropy varEntropy entropy] = getEntropyStats(im, labels, ...
                                                          centerInfo);
    
    if (~shouldPrint)
        return
    end

    
    avgIntensity = getAvgIntensity(centerInfo);
    avgVol = getAvgVol(centerInfo);
    
    [varIntensity varIntensitysv] = getVarIntensity(centerInfo,im,labels);
    varVolume = getVarVolume(centerInfo);
    
    surfaceArea = getSurfaceArea(labels, centerInfo);
    avgSurfaceArea = getAvgSurfaceArea(surfaceArea);
    varSurfaceArea = getVarSurfaceArea(surfaceArea);                                          

    [avgSpreadX avgSpreadY avgSpreadZ spreadX spreadY spreadZ] = ...
        getSpreads(labels, centerInfo);
    
    %get the neighbors
    numNeb = 4;
    neighbors = getNeighbors(centerInfo,numNeb);
    
    if (~isnan(tissues))
        tissueInfo = getTissueInfo(labels, tissues, centerInfo);
        
        [percentsvGM percentsvWM percentsvCSF] = ...
            getTissuePercentages(tissueInfo);
    end
    
        
    csvfilename = strcat(filename(1:end-4),'.csv');
    outFile = fopen(filename, 'w');
    outCSV = fopen(csvfilename, 'w');
    
    fprintf(outFile,'//Patient %d \n',id);
    fprintf(outFile, 'numsupervoxels(patientid%d, %d).\n', ...
            id, size(centerInfo,1));
    fprintf(outFile, 'averageintensity(patientid%d, %f).\n', id, avgIntensity);
    fprintf(outFile, 'averagevolume(patientid%d, %f).\n', id, avgVol);
    fprintf(outFile, 'averagesurfacearea(patientid%d, %f).\n', ...
            id, avgSurfaceArea);
    fprintf(outFile, 'averageentropy(patientid%d, %f).\n', id, avgEntropy);    
    fprintf(outFile, 'varintensity(patientid%d, %f).\n', id, varIntensity);
    fprintf(outFile, 'varvolume(patientid%d, %f).\n', id, varVolume);
    fprintf(outFile, 'varsurfacearea(patientid%d, %f).\n', id, varSurfaceArea);
    fprintf(outFile, 'varentropy(patientid%d, %f).\n', id, ...
            varEntropy);
    fprintf(outFile, 'avgxspread(patientid%d, %f).\n', id, ...
            avgSpreadX);
    fprintf(outFile, 'avgyspread(patientid%d, %f).\n', id, ...
            avgSpreadY);
    fprintf(outFile, 'avgzspread(patientid%d, %f).\n', id, ...
            avgSpreadZ);
    
    
    if (~isnan(tissues))
        fprintf(outFile, ['percentagepredominatelygm(patientid%d, %f).\n'],...
                id, percentsvGM);
        fprintf(outFile, ['percentagepredominatelywm(patientid%d, %f).\n'],...
                id, percentsvWM);
        fprintf(outFile, ['percentagepredominatelycsf(patientid%d, %f).\n'],...
                id, percentsvCSF);
    end
    
    %fprintf('Printing graph of average intensities');
    %graphIntensities(centerInfo);
    
    for i = 1:size(centerInfo,1)
        fprintf(outFile,'x(patientid%d, sv%d, %f).\n', id, i, ...
                centerInfo(i,1));
        fprintf(outFile,'y(patientid%d, sv%d, %f).\n', id, i, ...
                centerInfo(i,2));
        fprintf(outFile,'z(patientid%d, sv%d, %f).\n', id, i, ...
                centerInfo(i,3));
        fprintf(outFile,'xspread(patientid%d, sv%d, %f).\n', id, i, ...
                spreadX(i));
        fprintf(outFile,'yspread(patientid%d, sv%d, %f).\n', id, i, ...
                spreadY(i));
        fprintf(outFile,'zspread(patientid%d, sv%d, %f).\n', id, i, ...
                spreadZ(i));
        if (~isnan(tissues))
            fprintf(outFile,'percentgm(patientid%d, sv%d, %f).\n', id, i, ...
                    tissueInfo(i,3)/centerInfo(i,5));
            fprintf(outFile,'percentwm(patientid%d, sv%d, %f).\n', id, i, ...
                    tissueInfo(i,4)/centerInfo(i,5));
            fprintf(outFile,'percentcsf(patientid%d, sv%d, %f).\n', id, i, ...
                    tissueInfo(i,2)/centerInfo(i,5));
        end
        fprintf(outFile,'avgintensitypersv(patientid%d, sv%d, %f).\n', ...
                id, i, centerInfo(i,4));
        fprintf(outFile,'volumepersv(patientid%d, sv%d, %f).\n', ...
                id, i, centerInfo(i,5));
        fprintf(outFile,'varintensitypersv(patientid%d, sv%d, %f).\n', ...
                id, i, varIntensitysv(i));
        fprintf(outFile,'surfaceareapersv(patientid%d, sv%d, %f).\n', ...
                id, i, surfaceArea(i));
        fprintf(outFile,'entropypersv(patientid%d, sv%d, %f).\n', ...
                id, i, entropy(i));
        fprintf(outCSV,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f', ...
                centerInfo(i,1),centerInfo(i,2),centerInfo(i,3), ...
                spreadX(i),spreadY(i),spreadZ(i), ...
                tissueInfo(i,3)/centerInfo(i,5), ...
                tissueInfo(i,4)/centerInfo(i,5), ...
                tissueInfo(i,2)/centerInfo(i,5),centerInfo(i,4), ...
                centerInfo(i,5), varIntensitysv(i), ...
                surfaceArea(i),entropy(i));
        for neb = 1:numNeb
            printStr = strcat('neighbor(patientid%d, sv%d, sv%d).\n');
            fprintf(outFile,printStr, id, i, neighbors(i, neb));
            fprintf(outCSV,',%f',neighbors(i,neb));
        end
        fprintf(outCSV,'\n');
    end

    fclose(outFile);
    fclose(outCSV);
end

function issv = isSurfaceVoxel(i, j, k, labels)
% isSurfaceVoxel - Return true if voxel has another voxel of a
% different label as a neighbor in any cartesian direction
    
    issv = false;
    if (i==1 || i==size(labels, 1) || j==1 || j==size(labels, 2) || ...
        k==1 || k==size(labels, 3))
        
        issv = true;
        return;
    end

    for sub_i = [i-1 i+1]
        if labels(sub_i,j,k) ~= labels(i,j,k)
            issv = true;
            return;
        end
    end
    for sub_j = [j-1 j+1]
        if labels(i,sub_j,k) ~= labels(i,j,k)
            issv = true;
            return;
        end
    end

    for sub_k = [k-1 k+1]
        if labels(i,j,sub_k) ~= labels(i,j,k)
            issv = true;
            return;
        end
    end
end

function graphIntensities(centerInfo)
    figure
    plot(centerInfo(:,4),'x');
    title('Average supervoxel intensities');
    ylabel('Intensity');
    xlabel('Supervoxel intensity');
end


function avgIntensity = getAvgIntensity(centerInfo)

    avgIntensity = mean(centerInfo(:, 4));    
end

function avgVol = getAvgVol(centerInfo)

    avgVol = mean(centerInfo(:, 5));
end

function [varIntensity varIntensitysv] = getVarIntensity(centerInfo,im,labels)

    varIntensitysv = zeros(size(centerInfo,1));
    
    for i = 1:size(centerInfo,1)
        varIntensitysv(i) = var(im(labels==i));
    end
    
    varIntensity = var(centerInfo(:, 4));
end

function varVolume = getVarVolume(centerInfo)

    varVolume = var(centerInfo(:, 5));
end


function surfaceArea = getSurfaceArea(labels, centerInfo)

    surfaceArea = zeros(size(centerInfo, 1), 1);
    for i= 1:size(labels, 1)
        for j = 1:size(labels, 2)
            for k = 1:size(labels, 3)
                
                if isSurfaceVoxel(i, j, k, labels)
                    surfaceArea(labels(i, j, k)) = surfaceArea(labels(i, ...
                                                                      j, k)) + 1;
                end
            end
        end
    end
end

function avgSurfaceArea = getAvgSurfaceArea(surfaceArea)

    avgSurfaceArea = mean(surfaceArea);

end

function varSurfaceArea = getVarSurfaceArea(surfaceArea)

    varSurfaceArea = var(surfaceArea);
end

function [percentsvGM percentsvWM percentsvCSF] = ...
        getTissuePercentages(tissueInfo)
    
    totGM = 0;
    totWM = 0;
    totCSF = 0;
    for i=1:size(tissueInfo, 1)
        
        switch tissueInfo(i, 5)
          case 2
            totCSF = totCSF + 1;
          case 3
            totGM = totGM + 1;
          case 4
            totWM = totWM + 1;
        end
    end
    
    percentsvGM = (totGM) / (totGM + totWM + totCSF);
    percentsvWM = (totWM) / (totGM + totWM + totCSF);
    percentsvCSF = (totCSF) / (totGM + totWM + totCSF);    
end

function [avgEntropy varEntropy entropy] = getEntropyStats(im, labels, centerInfo)
    
    entropy = zeros(size(centerInfo,1),1);
    
    for i = 1:size(centerInfo,1)
        numCenters = centerInfo(i,5);
        entropySpread = (hist(double(im(labels==i)),255)./ numCenters) ...
            .*log(hist(double(im(labels==i)),255)./numCenters);
        entropy(i) = -1 * sum(entropySpread(~isnan(entropySpread))); 
    end
    
    avgEntropy = mean(entropy);
    varEntropy = var(entropy);
end

function [avgSpreadX avgSpreadY avgSpreadZ spreadX spreadY spreadZ] = getSpreads(labels, centerInfo)
    
    spreads = zeros(size(centerInfo, 1), 3);
    
    for i=1:size(spreads, 1)
        linearLabels = find(labels == i);
        if (size(linearLabels, 1) == 0)
            continue;
        end
        
        [rows cols pages] = ind2sub(size(labels), linearLabels);
        
        spreads(i, 1) = max(rows) - min(rows);
        spreads(i, 2) = max(cols) - min(cols);
        spreads(i, 3) = max(pages) - min(pages);
    end
    
    spreadX = spreads(:,1);
    spreadY = spreads(:,2);
    spreadZ = spreads(:,3);

    avgSpreadX = mean(spreads(:, 1));
    avgSpreadY = mean(spreads(:, 2));
    avgSpreadZ = mean(spreads(:, 3));
end

function neighbors = getNeighbors(centerInfo,numNeb)
   
    neighbors = zeros(size(centerInfo,1),numNeb);
    centers = centerInfo(:, 1:3)';
    numSV = size(centerInfo,1);
    
    for i = 1:numSV
        dist = distEuclidean(repmat(centers(:, i), 1, numSV),centers);
        
        [s, O] = sort(dist, 'ascend');
        
        neighbors(i,:) = O(2:numNeb+1);
    end
end