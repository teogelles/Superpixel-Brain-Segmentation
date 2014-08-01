% MDMLDEMO Demonstrates how to run MDML (mirror descent for metric learning), 
% which learns a low-dimensional (sparse) Mahalanobis metric, given
% labeled pairs of points. A pair of points (x, z) is labeled 'similar' (y 
% = 1) or dissimilar (y = -1). Sparsity here is in terms of the eigen-
% spectrum of the metric matrix, M, which is used to measure the distance 
% betwee two points (x, z) as
%
% d(x, z) = sqrt( (x - z)' M (x - z) ) = sqrt(  (x - z)' V E V'(x - z) ),
% 
% i.e., diag(E) is sparse. Projection into this lower-dimensional space can
% be achieved via
%                 P(x) = V(:, I) * sqrt(E(I)) * x ,
% where the indices I indicate the dimensions corresponding to the largest
% q eigenvalues of M/E.
%
% The Original (2d+8d garbage) is 2 relevant + 8 garbage dimensions.  
%
% Gautam Kunapuli (gkunapuli@gmail.com)
% January 12, 2013
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function [] = ADNI_mdmlDemo(filename,groupname)

K = 3;    % Number of nearest neighbors to consider


%%Our files. Wanna run new files? Putem here

if ~exist('groupname','var')
    filename = ['/scratch/tgelles1/summer2014/slicExact125/features/' ...
                'AllPat.csv'];
    groupname = ['/scratch/tgelles1/summer2014/slicExact125/features/' ...
                 'AllPat_groups.csv'];
end

[rawData.Xtrn rawData.ytrn rawData.Xtst rawData.ytst] = loadAndSplit(filename,groupname,1,3);

% Construct tuples of supervised pairs of data, using K nearest neighbors                   
[Isim, Idis] = generateLabeledPairs(rawData.Xtrn, rawData.ytrn, K);
mdmlData.X = rawData.Xtrn;
mdmlData.I = [Isim; Idis];
mdmlData.y = [ones(size(Isim, 1), 1); -ones(size(Idis, 1), 1)];

% Randomize the order of data points
J = randperm(length(mdmlData.y));
mdmlData.I = mdmlData.I(J, :);
mdmlData.y = mdmlData.y(J);

% Do MDML using hinge + frobenius
settings = mdmlSettings('loss', 'hinge',...     % Loss function
                        'breg', 'frobenius',... % Bregman divergence
                        'eta', 0.25,...         % Learning rate
                        'rho', 0.05,...         % Trace-norm regularization
                        'k', K,...              % Num. nearest neighbors
                        'acc', 1e-6,...         % Math. Accuracy/tolerance
                        'verbose', 0,...        % Display-level verbosity
                        'debug', 0,...          % Perform sanity checks
                        'progbar', 1);          % Progress bar
                    
% Run metric learning 
mdmlResults = mdml(mdmlData, settings);

% Now project the training data into a new space, whose basis is given by 
% the two largest eigenvalues, for plotting purposes
Mhalf = mdmlResults.V(:, [10, 9]) * diag(sqrt(mdmlResults.E([10, 9])));
Xptrn = rawData.Xtrn * Mhalf; 

% Compute the kNN error of this approach 
ypred = kNNClassify(Mhalf, rawData.Xtrn, rawData.ytrn, rawData.Xtst, K);
err = 1 - mean(ypred' == rawData.ytst);        
fprintf('Test error on 2D data set (H+F) = %g%%.\n', err*100);

%%%%
% Do MDML using logistic + von Neumann
settings = mdmlSettings('loss', 'logistic', 'breg', 'vonneumann',...
                        'eta', 0.5, 'rho', 0.05, 'k', K, 'acc', 1e-6,...
                        'verbose', 0, 'progbar', 1);
                    
% Run metric learning 
mdmlResults = mdml(mdmlData, settings);

% Now project the training data and the test data into a new space,
% whose basis is given by the two largest eigenvalues
Mhalf = mdmlResults.V(:, [10, 9]) * diag(sqrt(mdmlResults.E([10, 9])));
Xptrn = rawData.Xtrn * Mhalf; 

% Compute the kNN error of this approach 
ypred = kNNClassify(Mhalf, rawData.Xtrn, rawData.ytrn, rawData.Xtst, K);
err = 1 - mean(ypred' == rawData.ytst);        
fprintf('Test error on 2D data set (L+VN) = %g%%.\n\n', err*100);


% %
% % 3d data set
% %
% % Generate the training and test data
% rawData = generate3ClassData(L);

% % Construct tuples of supervised pairs of data, using K nearest neighbors                   
% [Isim, Idis] = generateLabeledPairs(rawData.Xtrn, rawData.ytrn, K);
% mdmlData.X = rawData.Xtrn;
% mdmlData.I = [Isim; Idis];
% mdmlData.y = [ones(size(Isim, 1), 1); -ones(size(Idis, 1), 1)];

% % Randomize the order of data points
% J = randperm(length(mdmlData.y));
% mdmlData.I = mdmlData.I(J, :);
% mdmlData.y = mdmlData.y(J);

% %%%%
% % Do MDML using hinge + frobenius
% settings = mdmlSettings('loss', 'hinge', 'breg', 'frobenius',...
%                         'eta', 0.5, 'rho', 0.2, 'k', K, 'acc', 1e-6,...
%                         'verbose', 0, 'progbar', 1);
                    
% % Run metric learning 
% mdmlResults = mdml(mdmlData, settings);

% % Now project the training data and the test data into a new space,
% % whose basis is given by the two largest eigenvalues
% Mhalf = mdmlResults.V(:, [10, 9]) * diag(sqrt(mdmlResults.E([10, 9])));
% Xptrn = rawData.Xtrn * Mhalf; 
        
% % Plot the Original (2d+8d garbage), excluding the original 8 noisy dimensions
% plotDataSet(handles(1, 2), rawData.Xtrn(:, rawData.I), rawData.ytrn,...
%                                 'Original (2d+8d garbage)');
                    
% % Plot data projected onto space learned with MDML: hinge + frobenius                    
% plotDataSet(handles(2, 2), Xptrn, rawData.ytrn, 'Projected Data (F + H)');

% % Compute the kNN error of this approach 
% ypred = kNNClassify(Mhalf, rawData.Xtrn, rawData.ytrn, rawData.Xtst, K);
% err = 1 - mean(ypred == rawData.ytst);        
% fprintf('Test error on 3D data set (F + H) = %g%%.\n', err*100);

                    
% %%%%
% % Do MDML using logistic + von Neumann
% settings = mdmlSettings('loss', 'logistic', 'breg', 'vonneumann',...
%                         'eta', 0.5, 'rho', 0.05, 'k', K, 'acc', 1e-6,...
%                         'verbose', 0, 'progbar', 1);
                    
% % Run metric learning 
% mdmlResults = mdml(mdmlData, settings);

% % Now project the training data into a new space, whose basis is given by 
% % the two largest eigenvalues, for plotting purposes
% Mhalf = mdmlResults.V(:, [10, 9]) * diag(sqrt(mdmlResults.E([10, 9])));
% Xptrn = rawData.Xtrn * Mhalf; 

% % Plot data projected onto space learned with MDML: logistic + vonneumann 
% plotDataSet(handles(3, 2), Xptrn, rawData.ytrn, 'Projected Data (VN + L)');
                    
% % Compute the kNN error of this approach 
% ypred = kNNClassify(Mhalf, rawData.Xtrn, rawData.ytrn, rawData.Xtst, K);
% err = 1 - mean(ypred == rawData.ytst);        
% fprintf('Test error on 3D data set (L + VN) = %g%%.\n\n', err*100);

function [] = plotDataSet(h, X, y, t)
colors = {'r', 'k', 'b'};
% styles = {'o', 's', '^'};

for j = 1:3
    scatter(h, X(y == j, 1), X(y == j, 2), 28, 'o', colors{j}, 'filled');
    hold(h, 'on');
end
axis(h, 'square');
title(h, t, 'FontSize', 12);


function data = generate2ClassData(L)
% Generate a 2D-normal data set with the following parameters
M = [0, 0]';
C = [2, 0.2;
     0.2, 2];
 
Xdata = mvnrnd(M, C, 10*L);
D = sqrt(sum(Xdata.^2, 2));
I = cell(2, 1);

sum(D < 0.7)
sum(D > 1)
I{1} = find(D < 0.85); I{1} = I{1}(1:L/2);
I{2} = find(D > 1.25); I{2} = I{2}(1:L/2);
Xdata = [Xdata(I{1}, :); Xdata(I{2}, :)];

Ydata = ones(L, 1);
Ydata(L/2+1:L) = 2;

% Now add 8 additional garbage dimensions. A good metric learning approach
% that learns low-dimensional representations should be able to strip these
% during learning
X = randn(L, 10);           % Generate 10D garbage data
I = randperm(10);           % Pick two random dimensions from 1:10
I = sort(I(1:2), 'ascend');
X(:, I) = Xdata;           % Insert our data set into those dimensions
Xdata = X;

% Split into 70% training and 30% test sets
J = randperm(L);
Ltrn = 0.7*L;
Jtrn = J(1:Ltrn);      Jtst = J(Ltrn+1:L);
data.Xtrn = Xdata(Jtrn, :); data.Xtst = Xdata(Jtst, :);
data.ytrn = Ydata(Jtrn);    data.ytst = Ydata(Jtst);
data.I = I;

function data = generate3ClassData(L)
% Generate a 2D-normal data set with the following parameters
M = [0, 0]';
C = [2, 0.2;
     0.2, 1];
 
Xdata = mvnrnd(M, C, L);

% Get indices of three classes and save them 
I = cell(3, 1);
I{1} = find(Xdata(:, 1) <= 0.9 & Xdata(:, 1) >= -0.9 &...
            Xdata(:, 2) <= 0.9 & Xdata(:, 2) >= -0.9);

I{2} = setdiff(find(Xdata(:, 1) < 0), I{1});
I{3} = setdiff(find(Xdata(:, 1) > 0), I{1});

Ydata = ones(L, 1);
Ydata(I{2}) = 2;
Ydata(I{3}) = 3;

% Now add 8 additional garbage dimensions. A good metric learning approach
% that learns low-dimensional representations should be able to strip these
% during learning
X = randn(L, 10);           % Generate 10D garbage data
I = randperm(10);           % Pick two random dimensions from 1:10
I = sort(I(1:2), 'ascend');
X(:, I) = Xdata;           % Insert our data set into those dimensions
Xdata = X;

% Split into 70% training and 30% test sets
J = randperm(L);
Ltrn = 0.7*L;
Jtrn = J(1:Ltrn);      Jtst = J(Ltrn+1:L);
data.Xtrn = Xdata(Jtrn, :); data.Xtst = Xdata(Jtst, :);
data.ytrn = Ydata(Jtrn);    data.ytst = Ydata(Jtst);
data.I = I;