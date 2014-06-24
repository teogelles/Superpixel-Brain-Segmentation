function tissueInfo = getTissueInfo(labels, tissues, centerInfo)
    
    
    numCenters = size(centerTracker, 1);
    
    tissueInfo = zeros(numCenters, 5);
    
    
    for i=1:size(tissues, 1)
        for j=1:size(tissues, 2)
            for k=1:size(tissues, 3)
                
                tissueInfo(labels(i, j, k), tissues(i, j, k)) = ...
                    tissueInfo(labels(i, j, k), tissues(i, j, k)) ...
                    + 1;
            end
        end
    end
    
    
    for i=1:numCenters
        
        [maxTrash argmax] = max(tissueInfo(i, :));
        tissueInfo(i, 5) = argmax;
    end
end