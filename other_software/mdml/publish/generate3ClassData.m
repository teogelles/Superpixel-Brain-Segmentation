function data = generate3ClassData(L)
% Generate a 2D-normal data set with the following parameters
M = [0, 0]';
C = [2, 0.2;
     0.2, 1];
 
% We're going to generate data and then delete some of the generated data,
% so increase the data set slightly to account for these shenanigans
L = round(L*1.25); 
Xdata = mvnrnd(M, C, L);

% Get indices of three classes and save them 
I = cell(3, 1);

% Get everything in the square |x| <= 0.9 
I{1} = find(Xdata(:, 1) <= 0.9 & Xdata(:, 1) >= -0.9 &...
            Xdata(:, 2) <= 0.9 & Xdata(:, 2) >= -0.9);

% Get everything in the rectangle |x| <= [-1.25, 1.25, -1.1, 1.1]
I{2} = find(Xdata(:, 1) <= 1.25 & Xdata(:, 1) >= -1.25 &...
            Xdata(:, 2) <= 1.1 & Xdata(:, 2) >= -1.1);
        
% Strip the gap between them        
Xdata(setdiff(I{2}, I{1}), :) = [];
L = size(Xdata, 1);
        
% Now make three classes
I{1} = find(Xdata(:, 1) <= 0.9 & Xdata(:, 1) >= -0.9 &...
            Xdata(:, 2) <= 0.9 & Xdata(:, 2) >= -0.9);

I{2} = setdiff(find(Xdata(:, 1) < 0), I{1});
I{3} = setdiff(find(Xdata(:, 1) > 0), I{1});

Ydata = zeros(L, 1);
Ydata(I{1}) = 1;
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
Ltrn = round(0.7*L);
Jtrn = J(1:Ltrn);      Jtst = J(Ltrn+1:L);
data.Xtrn = Xdata(Jtrn, :); data.Xtst = Xdata(Jtst, :);
data.ytrn = Ydata(Jtrn);    data.ytst = Ydata(Jtst);
data.I = I;