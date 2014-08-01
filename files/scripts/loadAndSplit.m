function [train trainID test testID] = loadAndSplit(filename, ...
                                                    groupname,fold,numFolds)
% Simple program that loads in CSV files and splits them into a
% training and testing set based upon the fold and numfolds options
% filename = csv file with nxd matrix, where n is data points and d
%   is dimensions
% groupname = csv file with a nx1 matrix where all values are
%   labels of the data points' group
% fold = number of the fold, must be <= numFolds
% numFolds = number of total folds that will be done    
    
    totalVectors = csvread(filename);
    totalIDs = csvread(groupname);
    
    % Make sure that the data is sorted so that we can easily
    % choose the right proportion of data we can take from each of
    % the seperate groups
    if ~issorted(totalIDs)
        [totalIDs origOrder] = sort(totalIDs);
        totalVectors = totalVectors(origOrder,:);
    end
        
    [train trainID test testID] = ...
        getTrainingAndTesting(totalVectors,totalIDs,fold,numFolds);
    
end

function [train trainID test testID] = ...
        getTrainingAndTesting(totalVectors,totalIDs,fold,numFolds)
    differentGroups = unique(totalIDs);
    numGroups = length(differentGroups);
    
    numOfEach = zeros(numGroups,1);
    numTotal = size(totalVectors,1);
    
    %upper and lower fold ratios
    lfr = (fold-1)/numFolds;
    ufr = fold/numFolds;
    
    testInd = zeros(numTotal,1);
    
    for i = 1:numGroups
        totalSoFar = sum(numOfEach);
        numOfEach(i) = length(find(totalIDs == differentGroups(i)));
        testInd(floor(totalSoFar + lfr*numOfEach(i) + 1):...
                floor(totalSoFar + ufr*numOfEach(i))) = 1;
    end
    
    %make training index as not testing
    trainInd = ~testInd;
    %convert the booleans to indices
    trainInd = find(trainInd);
    testInd = find(testInd);

    % get testing
    test = totalVectors(testInd,:);
    testID = totalIDs(testInd);
    
    % get training
    train = totalVectors(trainInd,:);
    trainID = totalIDs(trainInd);
end 

    
    
    