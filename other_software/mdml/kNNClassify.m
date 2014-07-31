% KNNCLASSIFY  Performs (possibly multi-class) classification on test data
% given labeled training data, and the number of nearest neighbors to
% consider using the specified Mahalanobis distance.
%
% YTEST = kNNClassify(L, Xtrain, ytrain, Xtest, K)
%    L                  Mahalanobis projection matrix
%    Xtrain, ytrain     Labeled training data
%    ytest              Data to be classified
%    K                  Number of nearest neighbors to consider
%
%  version 3.7
%  Gautam Kunapuli (gkunapuli@gmail.com)
%  January 17, 2012
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function [ytest, confidence] = kNNClassify(L, Xtrain, ytrain, Xtest, K)

% Compute the Mahalanobis distance between the train and test points
D = mahalanobisDistance(L, Xtrain, Xtest);

% Sort the distances
[weights, neighbors] = sort(D);

% Keep only the indices of the K nearest neighbors (the first one is the
% data point itself and should be skipped)
neighbors = neighbors(1:K, :);
weights = 1 ./ weights(1:K, :);

% Determine the class of each column using weights
neighborLabels = ytrain(neighbors);

nTest = size(Xtest, 1);
classes = unique(ytrain);
nClasses = length(classes);

if K == 1
    % If there is only one neighbor
    ytest = neighborLabels;
    confidence = ones(numel(neighborLabels), 1);
elseif nClasses == 1
    % If there is only one training data class
    ytest = classes(ones(nTest, 1));
    confidence = ones(nTest, 1);
else
    % All other cases
    classWeights = zeros(nClasses, nTest);

    for i = 1:nTest
        uniqueLabels = unique(neighborLabels(:, i));
        for j = uniqueLabels'
            classWeights(classes == j, i) = sum(weights(neighborLabels(:, i) == j, i));
        end
    end

    % Get the class that has the largest sum of weights contributed to by
    % the nearest neighbors, as well as the sum of weights itself, which
    % represents a confidence
    [confidence, Itest] = max(classWeights);
    ytest = classes(Itest);
    ytest = ytest';
    
    % Normalize the confidence to be a number between 0 and 1 for each
    % particular data point
    confidence = confidence./ sum(classWeights);
end