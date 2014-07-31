% MAHALANOBISDISTANCE  Computes the Mahalanobis distance between points in
% the matrix X and the points in the matrix Z.
%
% D = mahalanobisDistance(L, X, Z)
%   L   N x M matrix, where M is rank of the projection space of the metric
%   X   Lx x N matrix, where Lx is the number of points in the matrix X 
%   Z   Lz x N matrix, where Lz is the number of points in the matrix Z
%   D   Lx x Lz matrix, whose (i, j)-th element is mahalanobis(X(i), Z(j))
%
%  version 2.0
%  Gautam Kunapuli (gkunapuli@gmail.com)
%  January 17, 2012
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function D = mahalanobisDistance(L, X, Z)

% Get some dimensions
K = size(L, 2);
[Mx, Nx] = size(X);
[Mz, Nz] = size(Z);
if Nx ~= Nz
    error('Matrices X and Z should have the same number of columns!');
end

% Project the data on to L
LX = X*L;
LZ = Z*L;

% Compute the pairwise distances
Ux = ones(Mx, K); % Unit matrix to redimension LZ
Uz = ones(Mz, K); % Unit matrix to redimension LX
D = sqrt( (LX.^2)*Uz' + Ux*(LZ.^2)' - 2*LX*LZ' );