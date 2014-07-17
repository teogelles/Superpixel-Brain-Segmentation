function findBraingroupPercentages()
    filebase =  ['/scratch/tgelles1/summer2014/ADNI_features/' ...
                 'CSV_NORM/'];
    groupfile = strcat(filebase,'small_groups.csv');
    resultsfile = strcat(filebase,'organized_small_results.csv');
    
    groups = csvread(groupfile);
    results = csvread(resultsfile);
    
    ADind = find(groups == 0);
    MCIind = find(groups == 1);
    CNind = find(groups == 2);
    
    ADresults = results(ADind,:);
    MCIresults = results(MCIind,:);
    CNresults = results(CNind,:);
    
    numClusters = max(results(:,1));
    
    clustersInAD = zeros(numClusters,1);
    clustersInMCI = zeros(numClusters,1);
    clustersInCN = zeros(numClusters,1);
    
    for cluster = 1:numClusters
        clustersInAD(cluster) = length(find(ADresults(:,1) == ...
                                            cluster));
        clustersInMCI(cluster) = length(find(MCIresults(:,1) == ...
                                            cluster));
        clustersInCN(cluster) = length(find(CNresults(:,1) == ...
                                            cluster));
    end
    
    totalADclusters = sum(clustersInAD(:));
    totalMCIclusters = sum(clustersInMCI(:));
    totalCNclusters = sum(clustersInCN(:));
    
    clustersInAD = clustersInAD/totalADclusters;
    clustersInMCI = clustersInMCI/totalMCIclusters;
    clustersInCN = clustersInCN/totalCNclusters;
    
    fprintf('Percentage clusters of each group in AD:\n');
    disp(clustersInAD);
    fprintf('Percentage clusters of each group in MCI:\n');
    disp(clustersInMCI);
    fprintf('Percentage clusters of each group in CN:\n');
    disp(clustersInCN);
end
