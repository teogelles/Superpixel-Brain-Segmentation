% EIGCUMLATIVE   Compute the number of eigenvalues required to achieve a
% cumulative energy specified by the threshold. For example, for a list of
% M eigenvalues E, if the sum of the first N reaches the threshold of 0.9 
% i.e., 90% of the total energy then EIGCUMULATIVE returns [N, sum(E(1:N)].
% This can be used to compute the first N principal components, for
% instance, such that these N components capture 90% of the total variance.
%
%  version 1.0
%  Gautam Kunapuli (gkunapuli@gmail.com)
%  April 4, 2012
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function [N, variance] = eigCumulative(E, threshold)

E = sort(E, 'descend');
traceE = sum(abs(E));
L = length(E);

for i = 1:L
    variance = sum(abs(E(1:i)));
    if (variance / traceE) >= threshold
        break;
    end
end

N = i;
variance = variance / traceE;