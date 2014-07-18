% We've changed this file to screw around with the segmentation

function RunSegmentation(fileName, numClusters, numNeighbors, ...
                         clusterType, graphType)
    
    FileDir = '/scratch/tgelles1/summer2014/ADNI_features/CSV_NORM/';

    Data = csvread(strcat(FileDir, fileName, '.csv'));
    [m, n, d] = size(Data);

    if d >= 2
        Data = (squeeze(Data))';
    end

    % convert to double and normalize to [0,1]
    Data = double(Data);
    Data = Data';
    Data = normalizeData(Data);

    % now for the clustering
    fprintf('Creating Similarity Graph...\n');
    SimGraph = SimGraph_NearestNeighbors(Data, numNeighbors, 1);

    try
        comps = graphconncomp(SimGraph, 'Directed', false);
        fprintf('- %d connected components found\n', comps);
    end

    fprintf('Clustering Data...\n');
    C = SpectralClustering(SimGraph, numClusters, clusterType);

    
    
    
    % convert and restore full size
    D = convertClusterVector(C);

    csvwrite(strcat(FileDir, fileName, '_out.csv'), [D Data']);

    
    
    
    % if (shouldPlot)
    %     % reshape indicator vector into m-by-n
    %     S = reshape(D, m, n);

    %     % choose colormap
    %     if k == 2
    %         map = [0 0 0; 1 1 1];
    %     else
    %         map = zeros(3, k);
            
    %         for ii = 1:k
    %             ind = find(D == ii, 1);
    %             map(:, ii) = rData(:, ind);
    %         end
            
    %         map = map';
    %     end

    %     % plot image
    %     set(gca, 'Position', [0 0 1 1], 'Units', 'Normalized');

    %     if isequal(markEdges, 1)
    %         imshow(Img, 'Border', 'tight');
            
    %         lS = label2rgb(S);
    %         BW = im2bw(lS, graythresh(lS));
    %         [B, L] = bwboundaries(BW, 'holes');

    %         hold on;
            
    %         for k = 1:length(B)
    %             boundary = B{k};
    %             plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
    %         end
            
    %         hold off;
    %     else
    %         imshow(S, map, 'Border', 'tight');
    %     end

    %     hold on;

    %     axis off;
    %     truesize;
    %     hold off;
    % end

    % clear all;
end