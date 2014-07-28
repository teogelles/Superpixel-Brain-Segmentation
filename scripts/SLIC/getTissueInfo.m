%% tissueInfo
% Authors: Teo Gelles & Andrew Gilchrist-Scott
%
% This file contians the getTissueInfo function, which obtains
% statistics related to the distribution of tissues in the brain as
% well as the SLIC supervoxelation.
function tissueInfo = getTissueInfo(labels, tissues, centerInfo)
% tissueInfo - Get info about tissues and SLIC supervoxelation
%
% @param labels - Matrix of SLIC labels of an image
% @param tissues - Matrix of tissue segmentation of an image
% @param centerInfo - Set of data for the centers of the SLIC
% supervoxelation
%
% @returns tissueInfo - Matrix of information
    
    numCenters = size(centerInfo, 1);
    
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