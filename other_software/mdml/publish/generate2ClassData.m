function data = generate2ClassData(L)
% Generate a 2D-normal data set with the following parameters
M = [0, 0]';
C = [2, 0.2;
     0.2, 2];
 
Xdata = mvnrnd(M, C, 10*L);
D = sqrt(sum(Xdata.^2, 2));
I = cell(2, 1);

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