%This is an attempt to combine the cluster percentages that we
%learned from Spectral clustering to create identifiable markers
%for whole brains

function [brainVectors brainIDs] = makeMRIvectors(fileName)
       
    if ~exist('fileSeparator','var')
        fileSeparator = 'small';
    end
    
    if ~exist('fold','var')
        fold = 1;
    end
    
    if ~exist('numFolds','var')
        numFolds = 5;
    end
    
    numSV = 120;
    
    [path filestub ext] = fileparts(fileName);
    
    %Find appropriate files
    resultsfile = [path '/' filestub '_clustered' ext];
    
    %convenient hack; remove
    fileSeparator = 'med';
        
    filebase = ['/scratch/tgelles1/summer2014/slicExact120/' ...
                'features/CSV_NORM/'];
    groupfile = strcat(filebase,fileSeparator,'_groups.csv');
    groups = csvread(groupfile);
    results = csvread(resultsfile);
    
    numBrains = size(results,1)/120;
    numClusters = max(results(:,1));
    
    [brainVectors brainIDs] = getBrainVectors(results,numBrains,numClusters, ...
                                                      numSV,groups);
    
end

function [brainVectors brainIDs] = getBrainVectors(results, ...
                                                   numBrains, ...
                                                   numClusters,numSV,groups)
    
    brainVectors = zeros(numBrains,numClusters);
    brainIDs = zeros(numBrains,1);
    brain_offset = 1;

    for brain = 1:numBrains
        clusterCount = zeros(numClusters,1);
        % We know that every brain has numSV supervoxels
        for SV = brain_offset:(brain_offset+numSV-1)
            clusterCount(results(SV,1)) = clusterCount(results(SV,1)) + 1;
        end
        % Assume all SV in the brain are properly labeled
        brainIDs(brain) = groups(SV);
        percents = clusterCount/numClusters;
        brainVectors(brain,:) = percents;
        brain_offset = brain_offset + numSV;
    end 
end

