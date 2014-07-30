disp('Loading iris data');
X = csvread('/scratch/tgelles1/summer2014/slicExact120/features/CSV_NORM/total_ADNI.csv');
y = csvread('/scratch/tgelles1/summer2014/slicExact120/features/CSV_NORM/tot_groups.csv');

disp('Running ITML');
num_folds = 2;
knn_neighbor_size = 4;
acc = CrossValidateKNN(y, X, @(y,X) MetricLearningAutotuneKnn(@ItmlAlg, y, X), num_folds, knn_neighbor_size);

disp(sprintf('kNN cross-validated accuracy = %f', acc));

