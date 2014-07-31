% We've changed this file to screw around with the segmentation,
% and now, it's our file! Muah ha ha ha ha!

function [C centers] = ADNIkmeans(FileName,k)

    global debug;
    dState = debug;
    debug = false; % Change this line to {dis,en}able debug statements
    % relates to any ddisp and dprintf statements
    
    if ~exist('FileName','var')
        FileName = ['/scratch/tgelles1/summer2014/slicExact120/' ...
                    'features/CSV_NORM/organized_med.csv'];
    end
    
    if ~strcmp(FileName(1),'/')
        FileName = strcat(['/scratch/tgelles1/summer2014/' ...
                           'slicExact120/features/CSV_NORM/'], FileName);
    end
        
    if ~exist('k','var')
        fprintf('Setting k to 40\n');
        k = 40; %number of clusters
    end
        
    saveData  = true;      % Whether or not to save the data once
                           % computed
    normalize = false;     % Whether or not to normalize our data
    
    % Prepare data
    Data = csvread(FileName);
    [m n d] = size(Data);
    Data = double(Data);
    if normalize
        Data = normalizeData(Data')';
    end
    
    if isequal(saveData, 1)
        [savePath, saveFile, ~] = fileparts(FileName);
        savePath = strcat(savePath,'/');
        if normalize
            csvwrite(strcat(savePath,saveFile,'_normed.nld'), Data);
        end
    end
    
    % actual kmeans portion
    [C centers] = kmeans(Data, k, 'start', 'cluster', ...
                 'EmptyAction', 'singleton');

    results = zeros(size(Data,1),size(Data,2) + 1);
    results(:,1) = C';
    results(:,2:end) = Data;
    
    if saveData
        csvwrite(strcat(savePath,saveFile, ...
                        '_klustered.csv'),results);
    end
        
    debug = dState;
end

