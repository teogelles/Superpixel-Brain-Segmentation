%This is an attempt to combine the cluster percentages that we
%learned from Spectral clustering to create identifiable markers
%for whole brains

function makeMRIvectors(fileSeparator)
   
    
    if ~exist('fileSeparator','var')
        fileSeparator = 'small';
    end
    
    filebase =  ['/scratch/tgelles1/summer2014/slicExact120/features/' ...
                 'CSV_NORM/'];
    groupfile = strcat(filebase,fileSeparator,'_groups.csv');
    resultsfile = strcat(filebase,'organized_',fileSeparator,'_clustered.csv');
    
    groups = csvread(groupfile);
    results = csvread(resultsfile);
    
    numBrains = size(results,1)/120;
    numClusters = max(resultsfile(:,1));
    
    brainVectors = getBrainVectors(results,numBrains,numClusters);
    
   
end

function brainVectors = getBrainVectors(results,numBrains,numClusters)
    
    brainVectors = zeros(numBrains,numClusters)
    brain_offset = 1;

    for brain = 1:numBrains
        clusterCount = zeros(numClusters,1);
        % We know that every brain has 120 SV
        for SV = brain_offset:(brain_offset+119)
            clusterCount(results(SV,1)) = clusterCount(results(SV,1)) + 1;
        end
        percents = clusterCount/numClusters;
        brainVectors(brain,:) = percents;
        brain_offset = brain_offset + 120;
    end 
end