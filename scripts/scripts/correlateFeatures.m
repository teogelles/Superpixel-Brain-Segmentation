% Simple program to see if any of our features are correlated

function correlateFeatures(featureFileName)
    
    if ~exist('featureFileName','var')
        featureFileName = ['/scratch/tgelles1/summer2014/' ...
                           'ADNI_features/CSV/total_ADNI_SV.csv'];
    end
    
    features = csvread(featureFileName);
    numFeat = size(features,2);
    
    [r,p] = corrcoef(features)  % Compute sample correlation and p-values.
    [i,j] = find(r>.8);  % Find significant correlations.
    %[i,j] = find(p<0.05);  % Find significant correlations.
    
    for feat = 1:numFeat
        featIndex = find(i == feat);
        fprintf('\nFeature %d correlates strongly with:\n',feat);
        if size(featIndex(:)) == 0
            fprintf('Nothing\n');
            continue
        end
        j(featIndex)
    end
    fprintf('\nThat is all\n');
end
