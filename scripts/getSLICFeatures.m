function featureList = getSLICFeatures(im, labels, tissues, centerInfo, cropOffset, filename, pid)
    
    fprintf('Getting SLIC Features for Patient %d\n', pid);
    
    outFile = fopen(filename, 'w');
    
    fprintf(outFile,'//Patient %d \n',1);
    
    avgIntensity = getAvgIntensity(centerInfo);
    avgVol = getAvgVol(centerInfo);
    
    [varIntensity varIntensitySV] = getVarIntensity(centerInfo,im,labels);
    varVolume = getVarVolume(centerInfo);
    
    surfaceArea = getSurfaceArea(labels, centerInfo);
    avgSurfaceArea = getAvgSurfaceArea(surfaceArea);
    varSurfaceArea = getVarSurfaceArea(surfaceArea);                                          
    [avgEntropy varEntropy entropy] = getEntropyStats(im, labels, centerInfo);
    [avgSpreadX avgSpreadY avgSpreadZ spreadX spreadY spreadZ] = ...
        getSpreads(labels, centerInfo);
    
    if (~isnan(tissues))
        tissueInfo = getTissueInfo(labels, tissues, centerInfo);

        [percentSVGM percentSVWM percentSVCSF] = ...
            getTissuePercentages(tissueInfo);
    end
    
    fprintf(outFile, 'numSupervoxels(patientid%d, %d).\n', ...
            pid, size(centerInfo,1));
    fprintf(outFile, 'AverageIntensity(patientid%d, %f).\n', pid, avgIntensity);
    fprintf(outFile, 'AverageVolume(patientid%d, %f).\n', pid, avgVol);
    fprintf(outFile, 'AverageSurfaceArea(patientid%d, %f).\n', ...
            pid, avgSurfaceArea);
    fprintf(outFile, 'AverageEntropy(patientid%d, %f).\n', pid, avgEntropy);    
    fprintf(outFile, 'VarIntensity(patientid%d, %f).\n', pid, varIntensity);
    fprintf(outFile, 'VarVolume(patientid%d, %f).\n', pid, varVolume);
    fprintf(outFile, 'VarSurfaceArea(patientid%d, %f).\n', pid, varSurfaceArea);
    fprintf(outFile, 'VarEntropy(patientid%d, %f).\n', pid, ...
            varEntropy);
    fprintf(outFile, 'AvgXSpread(patientid%d, %f).\n', pid, ...
            avgSpreadX);
    fprintf(outFile, 'AvgYSpread(patientid%d, %f).\n', pid, ...
            avgSpreadY);
    fprintf(outFile, 'AvgZSpread(patientid%d, %f).\n', pid, ...
            avgSpreadZ);
    
    
    if (~isnan(tissues))
        fprintf(outFile, ['PercentagePredominatelyGM(patientid%d, %f).\n'],...
                pid, percentSVGM);
        fprintf(outFile, ['PercentagePredominatelyWM(patientid%d, %f).\n'],...
                pid, percentSVWM);
        fprintf(outFile, ['PercentagePredominatelyCSF(patientid%d, %f).\n'],...
                pid, percentSVCSF);
    end
    %fprintf('Printing graph of average intensities');
    %graphIntensities(centerInfo);
    
    % For now we're just going to print info for the first ten
    % supervoxels as a proof of concept formatting
    for i = round(size(centerInfo,1)/2):round(size(centerInfo,1)/2) ...
        + 10%1:size(centerInfo,1)
        fprintf(outFile,'x(patientid%d, SV%d, %f).\n', pid, i, ...
                centerInfo(i,1));
        fprintf(outFile,'y(patientid%d, SV%d, %f).\n', pid, i, ...
                centerInfo(i,2));
        fprintf(outFile,'z(patientid%d, SV%d, %f).\n', pid, i, ...
                centerInfo(i,3));
        fprintf(outFile,'xSpread(patientid%d, SV%d, %f).\n', pid, i, ...
                spreadX(i));
        fprintf(outFile,'ySpread(patientid%d, SV%d, %f).\n', pid, i, ...
                spreadY(i));
        fprintf(outFile,'zSpread(patientid%d, SV%d, %f).\n', pid, i, ...
                spreadZ(i));
        fprintf(outFile,'percentGM(patientid%d, SV%d, %f).\n', pid, i, ...
                tissueInfo(i,3)/centerInfo(i,5));
        fprintf(outFile,'percentWM(patientid%d, SV%d, %f).\n', pid, i, ...
                tissueInfo(i,4)/centerInfo(i,5));
        fprintf(outFile,'percentCSF(patientid%d, SV%d, %f).\n', pid, i, ...
                tissueInfo(i,2)/centerInfo(i,5));
        fprintf(outFile,'AvgIntensityPerSV(patientid%d, SV%d, %f).\n', ...
                pid, i, centerInfo(i,4));
        fprintf(outFile,'VolumePerSV(patientid%d, SV%d, %f).\n', ...
                pid, i, centerInfo(i,5));
        fprintf(outFile,'VarIntensityPerSV(patientid%d, SV%d, %f).\n', ...
                pid, i, varIntensitySV(i));
        fprintf(outFile,'SurfaceAreaPerSV(patientid%d, SV%d, %f).\n', ...
                pid, i, surfaceArea(i));
        fprintf(outFile,'EntropyPerSV(patientid%d, SV%d, %f).\n', ...
                pid, i, entropy(i));
    end
    
    
    
    % Feature cell array that can be passed out of this function
    
    if (~isnan(tissues))
        featureList = cell(12, 2);
        startTissueIndex = 10;
    else
        featureList = cell(9, 2);
    end
    
    featureList{1, 1} = 'Average Intensity';
    featureList{1, 2} = avgIntensity;
    featureList{2, 1} = 'Average Volume';
    featureList{2, 2} = avgVol;
    featureList{3, 1} = 'Average Surface Area';
    featureList{3, 2} = avgSurfaceArea;
    featureList{4, 1} = 'Intensity Variance';
    featureList{4, 2} = varIntensity;
    featureList{5, 1} = 'Volume Variance';
    featureList{5, 2} = varVolume;
    featureList{6, 1} = 'Surface Area Variance';
    featureList{6, 2} = varSurfaceArea;
    featureList{7, 1} = 'Average Spread over x-axis';
    featureList{7, 2} = avgSpreadX;
    featureList{8, 1} = 'Average Spread over y-axis';
    featureList{8, 2} = avgSpreadY;
    featureList{9, 1} = 'Average Spread over z-axis';
    featureList{9, 2} = avgSpreadZ;
    
    if (~isnan(tissues))
        featureList{startTissueIndex, 1} = ['Percentage of Predominately GM ' ...
                            'Supervoxels'];
        featureList{startTissueIndex, 2} = percentSVGM;
        featureList{startTissueIndex+1, 1} = ['Percentage of Predominately WM ' ...
                            'Supervoxels'];
        featureList{startTissueIndex+1, 2} = percentSVWM;
        featureList{startTissueIndex+2, 1} = ['Percentage of Predominately CSF ' ...
                            'Supervoxels'];
        featureList{startTissueIndex+2, 2} = percentSVCSF;
    end
    
    fclose(outFile);
end

function isSV = isSurfaceVoxel(i, j, k, labels)
    
    
    isSV = false;
    if (i==1 || i==size(labels, 1) || j==1 || j==size(labels, 2) || ...
        k==1 || k==size(labels, 3))
        
        isSV = true;
        return;
    end

    for sub_i = [i-1 i+1]
        if labels(sub_i,j,k) ~= labels(i,j,k)
            isSV = true;
            return;
        end
    end
    for sub_j = [j-1 j+1]
        if labels(i,sub_j,k) ~= labels(i,j,k)
            isSV = true;
            return;
        end
    end

    for sub_k = [k-1 k+1]
        if labels(i,j,sub_k) ~= labels(i,j,k)
            isSV = true;
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

function [varIntensity varIntensitySV] = getVarIntensity(centerInfo,im,labels)

    varIntensitySV = zeros(size(centerInfo,1));
    
    for i = 1:size(centerInfo,1)
        varIntensitySV(i) = var(im(labels==i));
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

function [percentSVGM percentSVWM percentSVCSF] = ...
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
    
    percentSVGM = (totGM) / (totGM + totWM + totCSF);
    percentSVWM = (totWM) / (totGM + totWM + totCSF);
    percentSVCSF = (totCSF) / (totGM + totWM + totCSF);    
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