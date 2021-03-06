% The following code is using an SVM to attempt to classify our
% brains with a matrix that has all the bain info concatenated onto
% a single line, thus the brain block.
%
% Written by Teo Gelles and Andrew Gilchrist-Scott
%
% Note: this used libsvm rather than matlab's built-in SVM software

function accuracy = brainBlockSVM(fileName, groupName, numFolds)
    
    if ~exist('numFolds','var')
        numFolds = 2;
    end
    
    if ~exist('groupName','var')
        groupName = ['/scratch/tgelles1/summer2014/slicExact125/' ...
                    'features/EqAllPat_groups.csv'];
    end
    
    if ~exist('fileName','var')
        fileName = ['/scratch/tgelles1/summer2014/slicExact125/' ...
                    'features/EqAllPat.csv'];
    end
    
    brainVectors = csvread(fileName);
    
    % These lines make sure that the data is normalized such that
    % each column is on a range of 0 to 1
    brainVectors = (brainVectors - repmat(min(brainVectors,[],1), ...
                           size(brainVectors,1),1))*spdiags(1./ ...
                                                      (max(brainVectors,[],1)-min(brainVectors,[],1))',0,size(brainVectors,2),size(brainVectors,2));
    brainVectors = brainVectors(:,all(~isnan(brainVectors)));
    
    % disp(brainVectors(:,1:15))
    % pause;

    brainIDs = csvread(groupName);
    
    totalAcc = zeros(numFolds,1);
    for fold = 1:numFolds
        [train trainID test testID] = ...
            getTrainingAndTesting(brainVectors,brainIDs, fold,numFolds);
        
        % Trade off commented regions to try to vary parameters for
        % the SVM
        model = svmtrain(trainID, train);

        % bestcv = 0;
        % for log2c = -1:10,
        %     for log2g = -4:-4,
        %         cmd = ['-v 5 -c ', num2str(10^log2c), ' -g ', ...
        %                num2str(2^log2g), ' -t 2'];
        %         cv = svmtrain(trainID, train, cmd);
        %         if (cv > bestcv),
        %             bestcv = cv; bestc = 10^log2c; bestg = 2^log2g;
        %         end
        %         fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
        %     end
        % end

        % cmd = ['-c ', num2str(bestc), ' -g ', num2str(bestg), '-t 2'];
        % model = svmtrain(trainID, train, cmd);
        
        [predictions, acc, probs] = svmpredict(testID, test, ...
                                               model);
        
        makeConfusionMatrix(predictions,testID)
        
        %        disp([predictions, testID]);        
        totalAcc(fold) = acc(1);
    end
    
    fprintf('\n');

    accuracy = mean(totalAcc);
    maxaccuracy = max(totalAcc);
    minaccuracy = min(totalAcc);
    sdaccuracy = std(totalAcc);
    
    fprintf('Accuracy = %f\n',accuracy);
    fprintf('Max Accuracy = %f\n', maxaccuracy);
    fprintf('Min Accuracy = %f\n', minaccuracy);
    fprintf('SD Accuracy = %f\n', sdaccuracy);
end

function [train trainID test testID] = ...
        getTrainingAndTesting(brainVectors,brainIDs,fold,numFolds)
    
    numAD = length(find(brainIDs == 1));
    numMCI = length(find(brainIDs == 2));
    numCN = length(find(brainIDs == 3));
    
    fprintf('AD: %d, MCI: %d, CN: %d\n', numAD, numMCI, numCN);
    numBrains = numAD + numMCI + numCN;
    
    %upper and lower fold ratios
    lfr = (fold-1)/numFolds;
    ufr = fold/numFolds;
    
    testInd = zeros(numBrains,1);
    testInd(floor(lfr*numAD+1):floor(ufr*numAD)) = 1;
    testInd(floor(numAD + lfr*numMCI+1):floor(numAD + ufr*numMCI)) = 1;
    testInd(floor(numAD + numMCI + lfr*numCN+1):floor(numAD + numMCI + ...
                                                      ufr*numCN)) = ...
        1;
    
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

function makeConfusionMatrix(predictions,testID)
    k = length(unique(testID));
    conf = zeros(k);
    for i = 1:k
        for j = 1:k
            conf(i,j) = sum((predictions == i) .* (testID == j));
        end
    end
    fprintf('Confusion matrix:\n')
    disp(conf)
end

    
