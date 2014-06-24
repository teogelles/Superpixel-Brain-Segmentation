function featureList = getSLICFeatures(labels, centers, centerTracker, ...
                                               tissueFilename, ...
                                               indexList, filename, ...
                                               res)
    
    
    fprintf('Getting SLIC Features\n');
    
    avgIntensity = getAvgIntensity(centerTracker);
    avgVol = getAvgVol(centerTracker);
    
    varIntensity = getVarIntensity(centerTracker);
    varVolume = getVarVolume(centerTracker);
        
    surfaceArea = getSurfaceArea(centers, labels);
    avgSurfaceArea = getAvgSurfaceArea(surfaceArea);
    varSurfaceArea = getVarSurfaceArea(surfaceArea);

    tissueTracker = getTissueTracker(labels, centerTracker, ...
                                             tissueFilename, indexList, ...
                                             res);
    

    [percentSVGM percentSVWM percentSVCSF] = ...
        getTissuePercentages(tissueTracker);
        
    
    fprintf('Average Intensity: %f\n', avgIntensity);
    fprintf('Average Volume: %f\n', avgVol);
    fprintf('Average Surface Area: %f\n', avgSurfaceArea);
    fprintf('Intensity Variance: %f\n', varIntensity);
    fprintf('Volume Variance: %f\n', varVolume);
    fprintf('Surface Area Variance: %f\n', varSurfaceArea);
    fprintf(['Percentage of Predominately GM Supervoxels: ' ...
    '%f\n'], percentSVGM);
    fprintf(['Percentage of Predominately WM Supervoxels: ' ...
             '%f\n'], percentSVWM);
    fprintf(['Percentage of Predominately CSF Supervoxels: ' ...
             '%f\n'], percentSVCSF);
    %fprintf('Printing graph of average intensities');
    %graphIntensities(centerTracker);
    
    featureList = cell(6, 2);
    featureList{1, 1} = 'Average Intensity';
    featureList{1, 2} = avgIntensity;
    featureList{2, 1} = 'Avererage Volume';
    featureList{2, 2} = avgVol;
    featureList{3, 1} = 'Average Surface Area';
    featureList{3, 2} = avgSurfaceArea;
    featureList{4, 1} = 'Intensity Variance';
    featureList{4, 2} = varIntensity;
    featureList{5, 1} = 'Volume Variance';
    featureList{5, 2} = varVolume;
    featureList{6, 1} = 'Surface Area Variance';
    featureList{6, 2} = varSurfaceArea;

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

function graphIntensities(centerTracker)
    figure
    plot(centerTracker(:,4),'x');
    title('Average supervoxel intensities');
    ylabel('Intensity');
    xlabel('Supervoxel intensity');
end


function avgIntensity = getAvgIntensity(centerTracker)

    avgIntensity = mean(centerTracker(:, 4));    
end

function avgVol = getAvgVol(centerTracker)

    avgVol = mean(centerTracker(:, 5));
end
    
function varIntensity = getVarIntensity(centerTracker)

    varIntensity = var(centerTracker(:, 4));
end
 
function varVolume = getVarVolume(centerTracker)

    varVolume = var(centerTracker(:, 5));
end


function surfaceArea = getSurfaceArea(centers, labels)

    surfaceArea = zeros(size(centers, 1), 1);
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
        getTissuePercentages(tissueTracker)
    
    totGM = 0;
    totWM = 0;
    totCSF = 0;
    for i=1:size(tissueTracker, 1)
        
        switch tissueTracker(i, 5)
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