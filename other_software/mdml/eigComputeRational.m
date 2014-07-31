% EIGCOMPUTERATIONAL   Compute rationally, the eigenvalues of a rank-one 
% perturbation of a matrix with known eigenvalues. Given a matrix A with
% known eigenvalues E, EIGCOMPUTERATIONAL efficiently calculates the
% eigenvalues of the perturbed matrix, Ap = (A + rho*u*u'), where rho is 
% some scalar and u is a vector. 
% 
% NOTE: The algorithm takes t = V'*u as the input rather than z directly. 
%  i.e., we assume that Ap is expressed as Ap = V*(diag(E) + rho*t*t')*V'
%
% [mu, iters, execTime, status] = eigComputeRational(i, E, t, rho)
% returns a value mu such that the i-th eigenvalue of the pertubred matrix
% Ap can be computed as Ep(i) = E(i) + rho*mu. 
%
% Additonal parameters such as tolerance and maximum number of iterations
% can be specified: eigComputeRational(i, E, t, rho, tol, maxIter).
%
% The status values returned are
%   0 - No problems
%   1 - mu < 0 was encountered in some iteration
%   2 - mu > delta(i+1) was encountered in some iteration
%   3 - Maximum number of iterations exceeded 
%
% version 1.6
% Gautam Kunapuli (gkunapuli@gmail.com)
% April 5, 2012
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function [mu, iters, execTime, status] =...
                  eigComputeRational(i, E, t, rho, tol, maxIter)

if nargin < 6
    maxIter = 20;
end

if nargin < 5
    tol = 1e-10;
end

if nargin < 4
    error('Insufficient number of parameters provided for the method!');
end

% Start timing
tic;
status = 0;

% Compute delta
delta = (E - E(i)) / rho;
N = length(delta);

% Initialize some looping variables
converged = false;
iters = 0;

% The case for i = N is different, since many terms drop out of the
% interpolation expression making its computation easier
if i == N
    % Initialize the value of mu
    mu = 1;
    
    while ~converged && iters < maxIter
        iters = iters + 1;

        % Compute some helper values
        [psi, phi, Dpsi] = evaluateRationals(mu, delta, t, i, N);
        mu = mu + ( (1 + psi) / Dpsi ) * psi;
      
        % Guard the mu values
        if mu < 0
            mu = tol;
            % fprintf('WARNING: mu < 0: Fixed to mu = +tol. (Iter = %d), Nth eigenvalue.\n', iters);
            status = 1;
        end
        
        w = 1 + phi + psi;

        if abs(w) <= tol*N*(1 + abs(psi) + abs(phi))
            converged = true;
        end
    end
else
    % Initialize the value of mu
    I = [1:i-1, i+2:N];
    trest = t(I)' * ( t(I) ./ (delta(I) - delta(i+1)) );

    b = delta(i+1) + (t(i)^2 + t(i+1)^2) / (1 + trest);
    c = (delta(i+1) * t(i)^2) / (1 + trest);

    mu1 = b/2 - sqrt(b^2 - 4*c) / 2;
    mu2 = b/2 + sqrt(b^2 - 4*c) / 2;

    % Pick the smallest value for mu > 0
    if mu1 < tol && mu2 > tol
        mu = mu2;
    elseif mu1 > tol && mu2 < tol
        mu = mu1;
    else
        mu = min(abs(mu1), abs(mu2));
    end
    
    % Now iterate by constructing interpolating rationals
    while ~converged && iters < maxIter
        if isinf(mu) || isnan(mu)
            keyboard;
        end
        
        iters = iters + 1;

        % Compute some helper values
        D = delta(i+1) - mu;
        [psi, phi, Dpsi, Dphi] = evaluateRationals(mu, delta, t, i, N);

        % fprintf('%d %g %g %g %g %g\n', iters, mu, psi, phi, Dpsi, Dphi);
        
        c = 1 + phi - D*Dphi;

        a = (D*(1 + phi) + psi^2/Dpsi) / c + psi/Dpsi;

        w = 1 + phi + psi;
        b = (D*w*psi) / (Dpsi*c);

        % Update mu
        mu = mu + 2*b / (a + sqrt(a^2 - 4*b));
        
        % Guard the mu values
        if mu < 0
            mu = tol;
            % fprintf('WARNING: mu < 0. Fixed to mu = +tol.\n');
            status = 1;
        end
        
        if mu > delta(i+1)
            mu = delta(i+1) - tol;
            % fprintf('WARNING: mu > delta(i+1). Fixed to mu = delta(i+1)-tol.\n');
            status = 2;
        end

        % Check for convergence
        if abs(w) <= tol*N*(1 + abs(psi) + abs(phi))
            converged = true;
        end
    end
end
execTime = toc;

if iters == maxIter
    status = 3;
end

function [psi, phi, Dpsi, Dphi] = evaluateRationals(x, delta, u, i, N)
iL = 1:i;
iU = (i+1):N;

psi = u(iL)' * ( u(iL) ./ (delta(iL) - x) );
phi = u(iU)' * ( u(iU) ./ (delta(iU) - x) );

Dpsi = u(iL)' * ( u(iL) ./ (delta(iL) - x).^2 );
Dphi = u(iU)' * ( u(iU) ./ (delta(iU) - x).^2 );
