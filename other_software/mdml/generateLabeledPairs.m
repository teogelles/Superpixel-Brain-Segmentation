% GENERATELABELEDDATAPAIRS  Generates tuples of supervised pairs of data 
%  for metric learning based on K-nearest neighbors. 
%
%  version 2.0
%  Gautam Kunapuli (gkunapuli@gmail.com)
%  April 1, 2014
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function [Isim, Idis] = generateLabeledPairs(X, labels, K)

% Get dimensions
[L, N] = size(X);

% First construct the similarity matrix of squared Euclidean distances
U = ones(L, N);
S = X.^2*U' + U*(X.^2)' - 2*(X*X');

% Sort the similarities in descending order so that we can select the K 
% furthest data points with the same label for learning. We want to bring
% these points closer together when projected.
[~, neighbors] = sort(S, 'ascend');
neighbors(1, :) = [];
nLabels = labels(neighbors);

% Sort the points of the same label as the training point by distance
I = cell(L, 1);
for i = 1:L
    I(i) = { neighbors(nLabels(:, i) == labels(i), i); };
end

% Now go through this list and pick similar pairs while ensuring that there
% are no symmetric repeats
Isim = zeros(K*L, 2);
for i = 1:L
    J = I{i};
    Lj = length(J);
    
    if ~isempty(J)
        if Lj >= K
            J = J(1:K);
            Isim((i-1)*K+1:i*K, :) = [i*ones(K, 1), J];
        else
            Isim((i-1)*K+1:(i-1)*K+Lj, :) = [i*ones(Lj, 1), J];
        end
    end
           
    if ~isempty(J)
        for j = J'
            Jtemp = I{j};
            Jtemp(Jtemp == i) = [];
            I(j) = { Jtemp };
        end
    end
end
Isim(Isim(:, 1) == 0, :) = [];

% Sort the similarities in ascending order so that we can select the k
% nearest data points with a different label for learning. We want to send
% these points farther when projected.
[~, neighbors] = sort(S, 'ascend');
neighbors(1, :) = [];
nLabels = labels(neighbors);

I = cell(L, 1);
for i = 1:L
    I(i) = { neighbors(nLabels(:, i) ~= labels(i), i); };
end

% Now go through this list and pick similar pairs while ensuring that there
% are no symmetric repeats
Idis = zeros(K*L, 2);
for i = 1:L
    J = I{i};
    Lj = length(J);
    
    if ~isempty(J)
        if Lj >= K
            J = J(1:K);
            Idis((i-1)*K+1:i*K, :) = [i*ones(K, 1), J];
        else
            Idis((i-1)*K+1:(i-1)*K+Lj, :) = [i*ones(Lj, 1), J];
        end
    end
    
    if ~isempty(J)
        for j = J'
            Jtemp = I{j};
            Jtemp(Jtemp == i) = [];
            I(j) = { Jtemp };
        end
    end
end
Idis(Idis(:, 1) == 0, :) = [];

if size(Isim, 1) == 1
    Isim = Isim';
end

if size(Idis, 1) == 1
    Idis = Idis';
end