% The following program will be our attempt to pipe our spectral
% clustering (SC) results into a SVM to diagnosis AD vs CN vs MCI
%
% Written by Teo Gelles and Andrew Gilchrist-Scott
%
% Note: this used libsvm rather than matlab's built-in SVM software

function accuracy = SCpipeSVM(fileName, numFolds, k, neb, sigma)
    
    if ~exist('numFolds','var')
        numFolds = 3;
    end
    
    if ~exist('fileName','var')
        fileName = ['/scratch/tgelles1/summer2014/slicExact120/' ...
                    'features/CSV_NORM/organized_med.csv'];
    end
    
    SpectrallyCluster(fileName,k,neb,sigma);
    
    [brainVectors brainIDs] = makeMRIvectors(fileName);
    
    totalAcc = zeros(numFolds,1);
    for fold = 1:numFolds
        [train trainID test testID] = ...
            getTrainingAndTesting(brainVectors,brainIDs, fold,numFolds);
        
        model = svmtrain(trainID, train);
        
        [predictions, acc, probs] = svmpredict(testID, test, ...
                                               model);
        
        totalAcc(fold) = acc(1);
    end
    
    accuracy = mean(totalAcc);
    fprintf('Accuracy = %f\n',accuracy);
    
end

function [train trainID test testID] = ...
        getTrainingAndTesting(brainVectors,brainIDs,fold,numFolds)
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

    % get testing
    test = brainVectors(testInd,:);
    testID = brainIDs(testInd);
    
    % get training
    train = brainVectors(trainInd,:);
    trainID = brainIDs(trainInd);
end 
