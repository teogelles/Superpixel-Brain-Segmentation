% simple function to test sigma

function testSigma()
    
    filename = 'organized_med.csv';
    fileID = 'med';
    
    maxGap = 0;
    bestSig = Inf;
    
    for sig = 0:.1:3
        fprintf('Sigma of %2.1f\n',sig);
        SpectrallyCluster(filename,sig);
        percents = findBraingroupPercentages(fileID);
        % assume percents{1} has all clusters
        for i = 1:length(percents{1})
            if ((percents{1}(i) <= percents{2}(i)) && (percents{2}(i) ...
                                                       <= percents{3}(i))) ...
                    || ((percents{1}(i) >= percents{2}(i)) && ...
                        (percents{2}(i) >= percents{3}(i)))
                gap = abs(percents{1}(i) - percents{2}(i));
                if gap > maxGap
                    fprintf('New Best\n');
                    maxGap = gap;
                    bestSig = sig;
                end
            end
        end
        fprintf('\n');
    end
    
    if (bestSig ~= Inf)
        fprintf('Best sig: %2.1f with gap of %f\n',bestSig, maxGap);
    else
        fprintf('No gap found\n');
    end
end

                