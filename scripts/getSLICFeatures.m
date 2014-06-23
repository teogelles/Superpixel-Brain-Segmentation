function getSLICFeatures(labels, centers, centerTracker, filename)
    
    
    fprintf('Getting SLIC Features\n');
    
    avgIntensity = mean(centerTracker(:, 4));    
    avgVol = mean(centerTracker(:, 5));

    varIntensity = var(centerTracker(:, 4));
    varVolume = var(centerTracker(:, 5));
    
    
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
    
    avgSurfaceArea = mean(surfaceArea);
    varSurfaceArea = var(surfaceArea);
    
    fprintf('Average Intensity: %f\n', avgIntensity);
    fprintf('Average Volume: %f\n', avgVol);
    fprintf('Average Surface Area: %f\n', avgSurfaceArea);
    fprintf('Intensity Variance: %f\n', varIntensity);
    fprintf('Volume Variance: %f\n', varVolume);
    fprintf('Surface Area Variance: %f\n', varSurfaceArea);
    fprintf('Printing graph of average intensities');
    graphIntensities(centerTracker);
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
