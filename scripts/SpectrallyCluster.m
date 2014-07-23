% We've changed this file to screw around with the segmentation,
% and now, it's our file! Muah ha ha ha ha!

function SpectrallyCluster(FileName)

    global debug;
    dState = debug;
    debug = false; % Change this line to {dis,en}able debug statements
    % relates to any ddisp and dprintf statements
    
    if ~exist('FileName','var')
        FileName = ['/scratch/tgelles1/summer2014/ADNI_features/' ...
                    'CSV_NORM/organized_small.csv'];
    end
    
    if ~strcmp(FileName(1),'/')
        FileName = strcat(['/scratch/tgelles1/summer2014/' ...
                           'ADNI_features/CSV_NORM/'], FileName);
    end

    k         = 8;          % Number of Clusters
    Neighbors = 10;         % Number of Neighbors
    saveData  = false;      % Whether or not to save the data once computed
    
    Data = csvread(FileName);
    [m n d] = size(Data);
    Data = double(Data);
    Data = normalizeData(Data');
    
    

    if isequal(saveData, 1)
        [savePath, saveFile, ~] = fileparts(FileName);
        
        csvwrite(strcat(savePath,'/',saveFile,'_normed.nld'), Data);
    end

    % now for the clustering
    fprintf('Creating Similarity Graph...\n');
    SimGraph = SimGraph_NearestNeighbors(Data, Neighbors, 1);
    
    try
        comps = graphconncomp(SimGraph, 'Directed', false);
        fprintf('- %d connected components found\n', comps);
    end

    if saveData
        save(strcat(savePath,'/',saveFile,'_simgraph.mat'), SimGraph);
    end
    
    fprintf('Clustering Data...\n');
    [C, ~, ~, centers] = SpectralClustering(SimGraph, k, 2);
    
    ddisp('this is C')
    ddisp(C);
    ddisp('C is over')
    
    % convert and restore full size
    D = convertClusterVector(C);
    
    if debug; pause; end
    ddisp(D)
    ddisp(size(D))

    if saveData
        results = zeros(size(Data',1),size(Data',2) + 1);
        results(:,1) = D;
        results(:,2:end) = Data';
        csvwrite(strcat(savepath,'/',saveFile, ...
                        '_clustered.csv'),results);
        clear results
    end
        
    debug = dState;
end

