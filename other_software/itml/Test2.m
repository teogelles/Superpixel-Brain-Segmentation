disp('Loading iris data');
X = csvread('/scratch/tgelles1/summer2014/slicExact125/features/AllPat.csv');
y = csvread('/scratch/tgelles1/summer2014/slicExact125/features/AllPat_groups.csv');

X = X(1:20:end, :);
y = y(1:20:end, :);

disp('Running ITML');
num_folds = 2;
knn_neighbor_size = 4;
acc = CrossValidateKNN(y, X, @(y,X) MetricLearningAutotuneKnn(@ItmlAlg, ...
                                                  y, X), num_folds, ...
                       knn_neighbor_size);


disp(sprintf('kNN cross-validated accuracy = %f', acc));

