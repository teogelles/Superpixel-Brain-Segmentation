% The following code is using an SVM to attempt to classify our
% brains with a matrix that has all the bain info concatenated onto
% a single line, thus the brain block.
%
% Written by Teo Gelles and Andrew Gilchrist-Scott
%
% Note: this used libsvm rather than matlab's built-in SVM software

function accuracy = brainBlockSVM(fileName, groupName, numFolds)
    
    if ~exist('numFolds','var')
        numFolds = 5;
    end
    
    if ~exist('groupName','var')
        groupName = ['/scratch/tgelles1/summer2014/slicExact125/' ...
                    'features/AllPat_groups.csv'];
    end
    
    if ~exist('fileName','var')
        fileName = ['/scratch/tgelles1/summer2014/slicExact125/' ...
                    'features/AllPat.csv'];
    end
    
    brainVectors = csvread(fileName);
    brainIDs = csvread(groupName);
    
    totalAcc = zeros(numFolds,1);
    for fold = 1:numFolds
        [train trainID test testID] = ...
            getTrainingAndTesting(brainVectors,brainIDs, fold,numFolds);
        
        model = svmtrain(trainID, train);
        
        [predictions, acc, probs] = svmpredict(testID, test, ...
                                               model);
        
        totalAcc(fold) = acc(1);
        disp([predictions testID]);
    end
    
    accuracy = mean(totalAcc);
    fprintf('Accuracy = %f\n',accuracy);
    
end

function [train trainID test testID] = ...
        getTrainingAndTesting(brainVectors,brainIDs,fold,numFolds)
    numAD = length(find(brainIDs == 1));
    numMCI = length(find(brainIDs == 2));
    numCN = length(find(brainIDs == 3));
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

    % get testing
    test = brainVectors(testInd,:);
    testID = brainIDs(testInd);
    
    % get training
    train = brainVectors(trainInd,:);
    trainID = brainIDs(trainInd);
end 
