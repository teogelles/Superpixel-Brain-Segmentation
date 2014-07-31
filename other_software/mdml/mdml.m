%  MDML   Perform mirror descent for metric learning given labeled data
%  triplets. MDML learns a Mahalanobis matrix M = L*L', which is possibly
%  low rank, and positive semi-definite. The space spanned by L is
%  low-dimensional and captures the separation of the data best. 
%
%  The approach can handle multi-class data sets by appropriately setting
%  labels y = 1 (similar) or y = -1 (dissimilar). See GENERATELABELEDPAIRS
%  help for more information. 
%
%  MDML takes several parameters which are set using MDMLSETTINGS. Notable
%  among these are the learning rate (eta) and the trace-norm
%  regularization parameter (rho). One difference from the paper is that in
%  the paper, shrinkage thresholding was done with eta*rho, while here it
%  is done direclty with rho (and thus assuming that the scale eta was
%  absorbed). This was mainly to decouple the two parameters for better
%  control and cross-validation.
%
%  For full details, see:
%  "Mirror Descent for Metric Learning: A Unified Approach",
%  G. Kunapuli & J. W. Shavlik, 
%  Proc. 23rd European Conference on Machine Learning (Machine Learning and 
%  Knowledge Discovery in Databases, LNCS 7523, P. A. Flach, T. de Bie & 
%  N. Cristianini), pp. 859-874 (2012).
%
%  Gautam Kunapuli
%  May 13, 2014
%
% This program comes with ABSOLUTELY NO WARRANTY; See GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function results = mdml(data, settings)

% The progress bar is only turned on if the value of verbose is set to
% zero. Verbose always take precedence.
if settings.verbose == 1
    useProgressBar = 0;
else
    useProgressBar = settings.progbar;
end

% Start a progress bar
if useProgressBar
    progStr = sprintf('Metric Learning with %s + %s',...
                                        settings.breg, settings.loss);
    h = waitbar(0, 'Initializing...', 'name', progStr);
end
    
% Use default settings if none are provided
if nargin < 2
    settings = comidMetricSettings;
    if settings.verbose > 0
        fprintf('\nUsing default settings:');
        disp(settings);
    end
else
    if settings.verbose > 0
        fprintf('\nUsing the following settings:');
        disp(settings);
    end
end

% Check the data
if isfield(data, 'X')
    N = size(data.X, 2); % # dimensions of the data
else
    error('cmData.X does not exist!');
end

if isfield(data, 'I')
    T = size(data.I, 1); % # pairs of points to be used for learning
else
    error('cmData.I does not exist!');
end

if isfield(data, 'y')
    if length(data.y) ~= T
        error('Number of labels does not match number of data pairs!');
    end
else
    error('cmData.y does not exist!');
end

% Initialize the metric learning variables; if values are provided,
% hotstart metric learning with the variables
if isfield(settings.init, 'mu')
    mu = settings.init.mu;
else
    mu = 1;
end

if isfield(settings.init, 'V')
    V = settings.init.V;
else
    V = eye(N);
end

if isfield(settings.init, 'E')
    E = settings.init.E;
else
    E = ones(N, 1);
end

% Set function to compute (most of) the gradient of the loss function. This
% function computes a scalar a(m(t)), which is a function of the margin at
% iteration t, m(t). The gradient of the loss can be expressed as below,
%       dL / d M  =  a(m(t))*y(t)*(x(t) - z(t))*(x(t) - z(t))', and
%       dL / d mu = -a(m(t))*y(t)
% given a pair of points at iteration t: x(t) and z(t), with label y(t)
if strcmp(settings.loss, 'hinge')
    gradLoss = @aHingeLoss;
elseif strcmp(settings.loss, 'squared')
    gradLoss = @aSquaredLoss;
elseif strcmp(settings.loss, 'exponential')
    gradLoss = @aExponentialLoss;
elseif strcmp(settings.loss, 'logistic')
    gradLoss = @aLogisticLoss;
end

% Set functions for different Bregman functions, which are subtly different
% for each one. This is in order to ensure that the most numerically stable
% solution is calculated by the rank-one eigenvalue solver
if strcmp(settings.breg, 'frobenius')
    bregmanUpdate = @updateFrobenius;
elseif strcmp(settings.breg, 'vonneumann')
    bregmanUpdate = @updateVonNeumann;
elseif strcmp(settings.breg, 'burg')
    bregmanUpdate = @updateBurg;
end

% Start debugMode if requested by the user
debugMode = settings.debug;
verbosity = settings.verbose;

if verbosity >= 0
        fprintf('Performing metric learning with %s + %s\n',...
                                        settings.breg, settings.loss);
end

% Time the learning
learningTimer = tic;
P = settings.nPasses;
numUpdates = 0;

% Start the iteration
for p = 1:P
    
    if verbosity > 0
        fprintf('Pass %d of %d\n', p, P);
        fprintf('-----------+--------------+--------------+---------------+--------\n');
        fprintf('   Iter    |    margin    |      a       |       mu      | zeros\n');
        fprintf('-----------+--------------+--------------+---------------+--------\n');
    end
    
    for t = 1:T
        if useProgressBar == 1
            if mod((p-1)*T + t, 250) == 0
                % Update the waitbar
                waitbar(((p-1)*T + t)/T/P, h,...
                    sprintf('%4.1f%% done (now, Pass %d of %d)',...
                    ((p-1)*T + t)/T/P*100, p, P));
            end
        end
        
        % Get the indices of the current pair of points
        i = data.I(t, 1);
        j = data.I(t, 2);
        
        % Really should not have this if the pairs were constructed correctly.
        % However, just check and skip to be safe
        if i == j
            continue;
        end
        
        % Turn on debugging to see where the algorithm fails. The most common
        % conditions are checked below
        if debugMode == 1
            if sum(isnan(E)) > 0
                error('Have Nans in E at interation %d', t);
            end
            
            if sum(sum(isnan(V))) ~= 0 || sum(sum(V > 1/settings.acc)) ~= 0
                error('Have NaNs or infs in V at iteration %d. Potential problems\nin the eigenvalue solver!', t);
            end
            
            orthV = norm(V'*V - eye(N));
            if orthV > settings.acc
                error('Lost orthogonality of V in iteration %d (error = %g).', t, orthV);
            end
            
            if sum(E) < settings.acc
                warning('MDML:zeroEigs',...
                    'All eigenvalues became zero after iteration %d', t-1);
            end
        end
        
        % Compute hinge loss gradients using the full eigendecomposition, i.e.,
        % loss(Mt)
        u = (data.X(i, :) - data.X(j, :))';
        if sum(abs(u)) < settings.acc
            continue;
        end
        
        w = V'*u;
        m = data.y(t) * (mu - w'*(E .* w));
        a = data.y(t) * gradLoss(m);
        
        % There is nothing to update if a = 0. If |a| is large enough, perfom
        % the appropriate update based on the chosen Bregman divergence
        if abs(a) > settings.acc
            eta = settings.eta / sqrt(t);
            [V, E, mu] = bregmanUpdate(V, E, mu, w, a, eta, settings.rho,...
                settings.acc);
            numUpdates = numUpdates + 1;
        else
            a = 0;
        end
        
        if verbosity > 0
            numZeros = length(find(E < settings.acc));
            fprintf(' % 7d   | %+5.4e | %+5.4e | %+5.4e |  %5d\n',...
                t, m, a, mu, numZeros);
        end
    end
    
    % Finish debugging this pass
    if debugMode
        disp('-----------+--------------+--------------+---------------+--------');
    end
    fprintf('After pass %d, # pairs seen = %d, # updates = %d.\n', p, T*p, numUpdates);
end

% Destroy the progress bar
if useProgressBar, close(h); end

% Return the results
results.V = V;
results.E = E;
results.mu = mu;
results.time = toc(learningTimer);
results.T = T;
results.U = numUpdates;

numZeros = length(find(E < settings.acc));
fprintf('Time        = %g seconds\n', results.time);
fprintf('#Zeros in E = %d / %d\n', numZeros, N);

%
% Helper functions that compute gradients for various losses, and updates
% for various Bregman divergences.
%

function a = aHingeLoss(m)
if (m < 1), a = 1; else a = 0; end

function a = aSquaredLoss(m)
a = max(1-m, 0);

function a = aExponentialLoss(m)
a = exp(-m);

function a = aLogisticLoss(m)
a = exp(-m) / (1 + exp(-m));

function [V, E, mu] = updateFrobenius(V, E, mu, w, a, eta, rho, acc)
% Compute the eigendecomposition of the perturbation
[V, E] = eigRankOneUpdate(V, E, w, -eta*a, acc);

% Shrink the eigenvalues and project onto the positive semi-definite cone
E(E < rho) = 0;

% Update the margin
mu = max(mu + eta*a, 1);

function [V, E, mu] = updateVonNeumann(V, E, mu, w, a, eta, rho, acc)
% Compute the log(E), for the non-zero eigenvalues
NZ = find(abs(E) > acc);
E(NZ) = log(E(NZ));

% Compute the eigendecomposition of the perturbation using the non-zero
% eigen-decomposition only
%if abs(eta*a) > acc
[V(:, NZ), E(NZ)] = eigRankOneUpdate(V(:, NZ), E(NZ), w(NZ), -eta*a, acc);
%end

% Shrink the eigenvalues and project onto the positive semi-definite cone
E(E < rho) = 0;

% Compute the inverse of gradPsi(E), for the vonNeumann div
E(NZ) = exp(E(NZ));

% Update the margin
mu = max(exp(log(mu) + eta*a), 1);

function [V, E, mu] = updateBurg(V, E, mu, w, a, eta, rho, acc)
% Compute the pseudo-inverse of the current matrix, i.e, invert only
% the non-zero eigenvalues, leaving the zero eigenvalues as they are
rhoNew = -eta*a / (1 + eta*a*w'*(E .* w));

% Compute the eigendecomposition of the perturbation
[V, E] = eigRankOneUpdate(V, E, w, rhoNew, acc);
N = length(E);

% Compute the inverse of gradPsi(E), for the Burg div
Z = find(abs(E) < acc); % Indices of the zero eigenvalues
NZ = setdiff(1:N, Z);   % Indices of the non-zero eigenvalues

E(Z) = 0;
E(NZ) =  1 ./ E(NZ);

% Shrink the eigenvalues
E(E < rho) = 0;

% Compute the inverse of gradPsi(E), for the Burg div
Z = find(abs(E) < acc); % Indices of the zero eigenvalues
NZ = setdiff(1:N, Z);   % Indices of the non-zero eigenvalues

E(Z) = 0;
E(NZ) = 1 ./ E(NZ);

% Project back
E(E < 0) = 0;
E = sort(E);

% Update the margin
if 1/mu - eta*a == 0
    muInt = mu;
else
    muInt = 1 / ( 1/mu - eta*a );
end
mu = max(muInt, 1);