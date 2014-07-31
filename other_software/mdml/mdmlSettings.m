% MDMLSETTINGS   Set the parameters for mdml with the defaults being 
% specified in the parantheses:
%         'acc'     (1e-12)    Convergence tolerance
%         'loss'   ('hinge')   Loss function {'hinge', 'logistic'}
%         'breg' ('frobenius') Bregman function {'frobenius', 'vonNeumann'}
%         'eta'     (1.0)      Learning rate ( > 0 )
%         'rho'     (0.5)      Trace norm regularization parameter
%         'k'        (4)       Number of nearest neighbors to consider
%         'debug'    (0)       Flag to turn on error checking (slow)
%         'timing'   (1)       Flag to turn on timing stats
%         'verbose'  (1)       Flag to turn on printing (slow)
%         'progbar'  (0)       Flag to turn on waitbar (if verbose = 1,
%                                                         this is 0 always)
%
% version 2.0
% Gautam Kunapuli (gkunapuli@gmail.com)
% April 1, 2012
%
% This program comes with ABSOLUTELY NO WARRANTY; See the GNU General Public
% License for more details. This is free software, and you are welcome to 
% modify or redistribute it.

function settings = mdmlSettings(varargin)

% Allowed values for loss and Bregman functions
allowedLosses = [{'hinge'}, {'squared'}, {'logistic'}, {'exponential'}];
allowedBregs  = [{'frobenius'}, {'vonneumann'}, {'burg'}];

% Create the struct of default settings
settings = struct('acc',    1e-12,...       % Convergence tolerance
                  'loss',   'hinge',...     % Loss function
                  'breg',   'frobenius',... % Bregman function
                  'eta',    1.0,...         % Learning rate
                  'rho',    0.5,...         % Tr. norm regulariz. parameter
                  'k',      4,...           % Number of nearest neighbors
                  'debug',  0,...           % Flag to turn on error checking
                  'timing', 1,...           % Flag to turn on timing stats
                  'verbose', 1,...          % Flag to turn on printing
                  'progbar', 0,...          % Flag to turn on progress bar
                  'nPasses', 1,...          % Num. times to pass through data
                  'init', []);              % Initalization structure to 
                                            % hotstart MDML with V, E or mu
              

% Check input arguments
if mod(nargin, 2) == 1
    disp('WARNING: Have an odd number of arguments! Using defaults.');
    return;
end 
                
% if nargin == 0
%     disp(sprintf('Using default settings.'));
%     return;
% end

% Parse the arguments
numParamsToSet = nargin / 2;
for n = 1:numParamsToSet
    fieldName = cell2mat(varargin(2*n - 1));
    fieldValue = cell2mat(varargin(2*n));
    
    % Do some argument checking
    if ~ischar(fieldName)
        fprintf('ERROR: Parameter #%d''s name is not a string!\n', n);
        continue;
    end

    % String arguments are checked here
    if strcmp(fieldName, 'loss') || strcmp(fieldName, 'breg')
        if ~ischar(fieldValue)
            fprintf('ERROR: Parameter ''%s'' should be a string!\n', fieldName);
            continue;
        end
        
        if strcmp(fieldName, 'loss') && sum(strcmp(fieldValue, allowedLosses)) == 0
            fprintf('ERROR: ''%s'' is not an allowed loss function (see help).\n', fieldValue);
            continue;
        end
        
        if strcmp(fieldName, 'breg') && sum(strcmp(fieldValue, allowedBregs)) == 0
            fprintf('ERROR: ''%s'' is not an allowed Bregman function (see help).\n', fieldValue);
            continue;
        end
    % Numeric arguments are checked here        
    else
        if ~isnumeric(fieldValue)
            fprintf('ERROR: Parameter ''%s'' should be a number!\n', fieldName);
            continue;
        end
        
        if fieldValue < 0
            fprintf('ERROR: The value for ''%s'' is negative (%g)!\n', fieldName, fieldValue);
            continue;
        end
    end
            
    % If everything checks out, set the value
    if isfield(settings, fieldName)
        settings.(fieldName) = fieldValue; 
    end
end
