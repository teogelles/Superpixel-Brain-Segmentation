% A program which simply takes the centers matched and finds the
% nearest center or neighbor

function goodPercent = simpleMatch(fileName,fold,numFolds)

    if ~exist('fileName','var')
        fileID = '/scratch/tgelles1/summer2014/slicExact120/features/CSV_NORM/organized_med.csv';
    end
    
    if ~any(fold == 1:numFolds)
        fprintf(['Error: cannot use that fold (%d). Fold must be an ' ...
                 'int between 1 and %d\n'], fold, numFolds);
        return
    end
    
    [brainVectors brainIDs] = makeMRIvectors(fileName);
    
    [centers test testID] = getTrainingAndTesting(brainVectors, ...
                                                  brainIDs,fold, ...
                                                  numFolds);
    %convert labels of 0:2 to 1:3 for indexing
    testID = testID + 1;
    
    guessedID = guessLabels(centers, test);
    
    fprintf('Right on %d\n',length(find(testID == guessedID)));
    fprintf('Wrong on %d\n',length(find(testID ~= guessedID)));
    
    goodPercent = length(find(testID == guessedID))/length(testID);
    
end

function [centers test testID] = getTrainingAndTesting(brainVectors,brainIDs, ...
                                                                    fold,numFolds)
    numAD = length(find(brainIDs == 0));
    numMCI = length(find(brainIDs == 1));
    numCN = length(find(brainIDs == 2));
    numBrains = numAD + numMCI + numCN;
    
    %upper and lower fold ratios
    lfr = (fold-1)/numFolds;
    ufr = fold/numFolds;
    
    testInd = zeros(numBrains,1);
    testInd(floor(lfr*numAD+1):floor(ufr*numAD)) = 1;
    testInd(floor(numAD + lfr*numMCI+1):floor(numAD + ufr*numMCI)) = 1;
    testInd(floor(numAD + numMCI + lfr*numCN+1):floor(numAD + numMCI + ...
                                                      ufr*numCN)) = 1;
    %make training index as not testing
    trainInd = ~testInd;
    %convert the booleans to indices
    trainInd = find(trainInd);
    testInd = find(testInd);
    
    test = brainVectors(testInd,:);
    testID = brainIDs(testInd);
    
    % calculate the centers of the clusters
    trainBrains = brainVectors(trainInd,:);
    trainInd = brainIDs(trainInd);
    
    numTrain = size(trainBrains,1);
    
    %AD center, then MCI, then CN
    centers = zeros(3,size(trainBrains,2));
    for xyz = 1:3
        section = find(trainInd == (xyz - 1));
        numSection = length(section);
        centers(xyz,:) = sum(trainBrains(section,:))/numSection;
    end
end 

function guessedID = guessLabels(centers, test)
    
    numBrains = size(test,1);
    guessedID = zeros(numBrains,1);
    for i = 1:numBrains
        brain = test(i,:);
        % technically for this to really be distance we would need
        % to take the square root, but since we're only concerned
        % with relationships, (ie a^2 < b^2 implies a < b) we won't bother
        distTest = sum((repmat(brain,3,1) - centers).^2,2);
        [~, guessedID(i)] = min(distTest);
    end
end