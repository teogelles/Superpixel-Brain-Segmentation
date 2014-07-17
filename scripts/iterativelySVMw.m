% Simple test to see if we get a global maximum on one weight

function iterativelySVMw(testInd,doPlots)
    
    if ~exist('testInd','var')
        testInd = 1;
    end
    
    if ~exist('doPlots','var')
        doPlots = false;
    end
    
    w = ones(14,1);
    
    iStart = 0;
    iEnd = 100;
    iStep = 5;
    numI = floor((iEnd - iStart)/iStep) + 1;
    saveDir = strcat('/scratch/tgelles1/summer2014/SVM_temp/SVM_iterTest_w',...
                     num2str(testInd),'_',num2str(iStart),'-',...
                     num2str(iStep),'-',num2str(iEnd),'.mat');

    if exist(saveDir,'file')
        testResults = load(saveDir);
        testResults = testResults.testResults;
    else
        
        testResults = zeros(numI,3);
        test_i = 0;
        
        for i = iStart:iStep:iEnd
            w(testInd) = i;
            runResults = test_svmCluster(w);
            test_i = test_i + 1;
            testResults(test_i,:) = [i runResults(1) runResults(3)];
        end
        
        save(saveDir,'testResults');
    end
    
    % most accurate weight = maw
    [maxAcc mawInd] = max(testResults(:,2));
    maw = testResults(mawInd,1);
    % least inaccurate weight = liw
    [minInacc liwInd] = min(testResults(:,2));
    liw = testResults(liwInd,1);
    
    fprintf('\nMost accurate is %d at %3.1f\n',maw, maxAcc);
    fprintf('Least inaccurate is %d at %3.1f\n\n',liw, minInacc);
        
    if doPlots
        figure(testInd)
        plot(testResults(:,1),testResults(:,2),'b')
        hold on;
        plot(testResults(:,2),testResults(:,3),'r')
        tit = sprintf('Accuracy of SVM with variation in weight %d',testInd);
        title(tit);
        ylabel('Percent accuracy');
        xlabel('Value of weight %d',testInd);
        legend('Accuracte SV labels','Inaccurate SV label');
        hold off;
    end
end