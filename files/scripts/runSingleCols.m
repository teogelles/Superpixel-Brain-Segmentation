function runSingleCols() 
    
    direc = ['/scratch/tgelles1/summer2014/slicExact120/features/' ...
             'CSV_NORM/allButOneData/'];
    goodness = zeros(17,4);
    goodness(:,1) = 1:17;
    for numClusters = 3:6
        for i = 3:17
            fileN = [direc 'organized_med-' num2str(i) '.csv'];
            SpectrallyCluster(fileN,numClusters);
            goodness(i,numClusters-2) = simpleMatch(fileN,1,2);
        end
    end
    disp(goodness);
end