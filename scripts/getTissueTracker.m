function tissueTracker = getTissueTracker(labels, centerTracker, ...
                                                  tissueFilename, ...
                                                  indexList, resx)
    
    
    tissues = load_nii(tissueFilename);
    tissues = tissues.img(1:res:end,1:res:end,1:res:end);
    max(max(max(tissues)))
    size(tissues)
    tissues = tissues(indexList(1, 1):indexList(1, 2), indexList(2, ...
                                                      1):indexList(2, ...
                                                      2), indexList(3, ...
                                                      1):indexList(3, 2));
    size(tissues)
    size(labels)
    
    numCenters = size(centerTracker, 1);
    
    tissueTracker = zeros(numCenters, 5);
    
    
    for i=1:size(tissues, 1)
        for j=1:size(tissues, 2)
            for k=1:size(tissues, 3)
                
                tissueTracker(labels(i, j, k), tissues(i, j, k)) = ...
                    tissueTracker(labels(i, j, k), tissues(i, j, k)) ...
                    + 1;
            end
        end
    end
    
    
    for i=1:numCenters
        
        [maxTrash argmax] = max(tissueTracker(i, :));
        tissueTracker(i, 5) = argmax;
    end
end