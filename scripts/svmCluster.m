% Attempt at using the MATLAB standard SVM library to do clustering
% on our ADNI data

function svmCluster(numAD,numMCI,numCN,numADtest,numMCItest, numCNtest);
    
    if (~exist('numAD','var')) || (~exist('numMCI','var')) || ...
            (~exist('numCN','var')) 
        numAD = 20;
        numMCI = 20;
        numCN = 20;
    end
    
    if (~exist('numADtest','var')) || (~exist('numMCItest','var')) || ...
            (~exist('numCNtest','var')) 
        numADtest = 5;
        numMCItest = 5;
        numCNtest = 5;
    end

    global trainMCI;
    trainMCI = true;
    
    csvFileHead = '/scratch/tgelles1/summer2014/ADNI_features/CSV/';
    
    [training groups] = getTraining(csvFileHead, numAD, numMCI, ...
                                                 numCN);
    if trainMCI
        [AD_SVMstruct CN_SVMstruct MCI_SVMstruct] = trainSVM(training,groups);
    else
        [AD_SVMstruct CN_SVMstruct] = trainSVM(training,groups);
    end
    
    startAD = numAD + 1;
    endAD = numAD + numADtest;
    startMCI = numMCI + 1;
    endMCI = numMCI + numMCItest;
    startCN = numCN + 1;
    endCN = numCN + numCNtest;
    
    [ADtest MCItest CNtest] = getTesting(csvFileHead, startAD, endAD, ...
                                                      startMCI, ...
                                                      endMCI, startCN, ...
                                                      endCN);
    svmTest(AD_SVMstruct,ADtest,MCI_SVMstruct,MCItest,CN_SVMstruct,CNtest);
    fprintf('Done\n');
end

function [training groups] = getTraining(csvFileHead,numAD,numMCI,numCN)

    fprintf('Getting training images...\n');
    
    global trainMCI;
    
    AD = cell(numAD,1);
    MCI = cell(numMCI,1);
    CN = cell(numCN,1);
    
    numADSV = 0;
    numMCISV = 0;
    numCNSV = 0;
    
    for i = 1:numAD
        % hack, fix sometime; for some reason AD012 doesn't exist
        if i == 12
            continue
        end
        filename = strcat(csvFileHead,'AD',sprintf('%03d',i), ...
                          '.csv');
        AD{i} = csvread(filename);
        AD{i} = removeBackgroundSV(AD{i});
        numADSV = numADSV + size(AD{i},1);
    end
    
    for i = 1:numMCI
        filename = strcat(csvFileHead,'MCI',sprintf('%03d',i), ...
                          '.csv');
        MCI{i} = csvread(filename);
        MCI{i} = removeBackgroundSV(MCI{i});
        numMCISV = numMCISV + size(MCI{i},1);
    end
    
    for i = 1:numCN
        filename = strcat(csvFileHead,'CN',sprintf('%03d',i), ...
                          '.csv');
        CN{i} = csvread(filename);
        CN{i} = removeBackgroundSV(CN{i});
        numCNSV = numCNSV + size(CN{i},1);
    end
    
    % Assume AD{1} has all features
    training = zeros(numADSV + numMCISV + numCNSV, size(AD{1},2));
    groupAD = cell(numADSV + numMCISV + numCNSV,1);
    if trainMCI
        groupMCI = cell(numADSV + numMCISV + numCNSV,1);
    end
    groupCN = cell(numADSV + numMCISV + numCNSV,1);
    if trainMCI
        groups  = cell(3,1);
    else
        groups = cell(2,1);
    end
    %for indexing in the large training matrix
    tot_sv = 1;
    
    for i = 1:numAD
        sizeOfThis = size(AD{i},1);
        training(tot_sv:(tot_sv + sizeOfThis-1),:) = AD{i};
        [groupAD{tot_sv:(tot_sv + sizeOfThis-1)}] = deal('AD');
        if trainMCI
            [groupMCI{tot_sv:(tot_sv + sizeOfThis-1)}] = deal(['not ' ...
                                'MCI']);
        end
        [groupCN{tot_sv:(tot_sv + sizeOfThis-1)}] = deal('not CN');
        tot_sv = tot_sv + sizeOfThis;
    end
    for i = 1:numMCI
        sizeOfThis = size(MCI{i},1);
        training(tot_sv:(tot_sv + sizeOfThis-1),:) = MCI{i};
        [groupAD{tot_sv:(tot_sv + sizeOfThis-1)}] = deal('not AD');
        if trainMCI
            [groupMCI{tot_sv:(tot_sv + sizeOfThis-1)}] = ...
                deal('MCI');
        end
        [groupCN{tot_sv:(tot_sv + sizeOfThis-1)}] = deal('not CN');
        tot_sv = tot_sv + sizeOfThis;
    end
    for i = 1:numCN
        sizeOfThis = size(CN{i},1);
        training(tot_sv:(tot_sv + sizeOfThis-1),:) = CN{i};
        [groupAD{tot_sv:(tot_sv + sizeOfThis-1)}] = deal('not AD');
        if trainMCI
            [groupMCI{tot_sv:(tot_sv + sizeOfThis-1)}] = deal(['not ' ...
                                'MCI']);
        end
        [groupCN{tot_sv:(tot_sv + sizeOfThis-1)}] = deal('CN');
        tot_sv = tot_sv + sizeOfThis;
    end
    
    groups{1} = groupAD;
    if trainMCI
        groups{2} = groupMCI;
        groups{3} = groupCN;
    else
        groups{2} = groupCN;
    end
end

function [ADtest MCItest CNtest] = getTesting(csvFileHead, startAD, endAD, ...
                                                           startMCI, ...
                                                           endMCI, startCN, ...
                                                           endCN)

    fprintf('Getting testing image(s)...\n');
    
    numAD =  (endAD - startAD) + 1;
    numMCI = (endMCI - startMCI) + 1;
    numCN = (endCN - startCN) + 1;
    
    
    AD = cell(numAD,1);
    MCI = cell(numMCI,1);
    CN = cell(numCN,1);
    
    numADSV = 0;
    numMCISV = 0;
    numCNSV = 0;
    
    for i = startAD:endAD
        filename = strcat(csvFileHead,'AD',sprintf('%03d',i), ...
                          '.csv');
        AD{i - startAD + 1} = csvread(filename);
        AD{i - startAD + 1} = removeBackgroundSV(AD{i - startAD + 1});
        numADSV = numADSV + size(AD{i - startAD + 1},1);
    end
    
    for i = startMCI:endMCI
        filename = strcat(csvFileHead,'MCI',sprintf('%03d',i), ...
                          '.csv');
        MCI{i - startMCI + 1} = csvread(filename);
        MCI{i - startMCI + 1} = removeBackgroundSV(MCI{i - startMCI + 1});
        numMCISV = numMCISV + size(MCI{i - startMCI + 1},1);
    end
    
    for i = startCN:endCN
        filename = strcat(csvFileHead,'CN',sprintf('%03d',i), ...
                          '.csv');
        CN{i - startCN + 1} = csvread(filename);
        CN{i - startCN + 1} = removeBackgroundSV(CN{i - startCN + 1});
        numCNSV = numCNSV + size(CN{i - startCN + 1},1);
    end
    
    % Assume first image has all features
    ADtest = zeros(numADSV,size(AD{1},2));
    MCItest = zeros(numMCISV,size(MCI{1},2));
    CNtest = zeros(numCNSV,size(CN{1},2));
    
    %for indexing in the testing martices
    AD_sv = 1;
    MCI_sv = 1;
    CN_sv = 1;
    
    for i = 1:numAD
        sizeOfThis = size(AD{i},1);
        ADtest(AD_sv:(AD_sv + sizeOfThis-1),:) = AD{i};
        AD_sv = AD_sv + sizeOfThis;
    end
    for i = 1:numMCI
        sizeOfThis = size(MCI{i},1);
        MCItest(MCI_sv:(MCI_sv + sizeOfThis-1),:) = MCI{i};
        MCI_sv = MCI_sv + sizeOfThis;
    end 
    for i = 1:numCN
        sizeOfThis = size(CN{i},1);
        CNtest(CN_sv:(CN_sv + sizeOfThis-1),:) = CN{i};
        CN_sv = CN_sv + sizeOfThis;
    end
end

function svmTestWithMCI(AD_SVMstruct,ADtest,MCI_SVMstruct,MCItest, ...
                        CN_SVMstruct,CNtest)
    
    fprintf('Classifying testing images');
    
    ADonAD = svmclassify(AD_SVMstruct,ADtest);
    fprintf('.');
    ADonMCI = svmclassify(MCI_SVMstruct,ADtest);
    fprintf('.');
    ADonCN = svmclassify(CN_SVMstruct,ADtest);
    fprintf('.');
    MCIonAD = svmclassify(AD_SVMstruct,MCItest);
    fprintf('.');
    MCIonMCI = svmclassify(MCI_SVMstruct,MCItest);
    fprintf('.');
    MCIonCN = svmclassify(CN_SVMstruct,MCItest);
    fprintf('.');
    CNonAD = svmclassify(AD_SVMstruct,CNtest);
    fprintf('.');
    CNonMCI = svmclassify(MCI_SVMstruct,CNtest);
    fprintf('.');
    CNonCN = svmclassify(CN_SVMstruct,CNtest);
    fprintf('.\n');
    
    total = size(ADtest,1) + size(MCItest,1) + size(CNtest,1);
    
    
    allRight = 0;
    kindaRight = 0;
    kindaWrong = 0;
    allWrong = 0;
    
    for i = 1:size(ADtest,1)
        if strcmp(ADonAD{i},'AD')
            if strcmp(ADonMCI{i},'not MCI') && strcmp(ADonCN{i},'not CN')
                allRight = allRight + 1;
            elseif strcmp(ADonMCI{i},'not MCI') || strcmp(ADonCN{i},'not CN')
                kindaRight = kindaRight + 1;
            else
                kindaWrong = kindaWrong + 1;
            end
        else
            if strcmp(ADonMCI{i},'not MCI') && strcmp(ADonCN{i},'not CN')
                kindaRight = kindaRight + 1;
            elseif strcmp(ADonMCI{i},'not MCI') || strcmp(ADonCN{i},'not CN')
                kindaWrong = kindaWrong + 1;
            else
                allWrong = allWrong + 1;
            end
        end
    end
    
    for i = 1:size(MCItest,1)
        if strcmp(MCIonAD{i},'not AD')
            if strcmp(MCIonMCI{i},'MCI') && strcmp(MCIonCN{i},'not CN')
                allRight = allRight + 1;
            elseif strcmp(MCIonMCI{i},'MCI') || strcmp(MCIonCN{i},'not CN')
                kindaRight = kindaRight + 1;
            else
                kindaWrong = kindaWrong + 1;
            end
        else
            if strcmp(MCIonMCI{i},'MCI') && strcmp(MCIonCN{i},'not CN')
                kindaRight = kindaRight + 1;
            elseif strcmp(MCIonMCI{i},'MCI') || strcmp(MCIonCN{i},'not CN')
                kindaWrong = kindaWrong + 1;
            else
                allWrong = allWrong + 1;
            end
        end
    end
    
    for i = 1:size(CNtest,1)
        if strcmp(ADonAD{i},'not AD')
            if strcmp(CNonMCI{i},'not MCI') && strcmp(CNonCN{i},'CN')
                allRight = allRight + 1;
            elseif strcmp(CNonMCI{i},'not MCI') || strcmp(CNonCN{i},'CN')
                kindaRight = kindaRight + 1;
            else
                kindaWrong = kindaWrong + 1;
            end
        else
            if strcmp(CNonMCI{i},'not MCI') && strcmp(CNonCN{i},'CN')
                kindaRight = kindaRight + 1;
            elseif strcmp(CNonMCI{i},'not MCI') || strcmp(CNonCN{i},'CN')
                kindaWrong = kindaWrong + 1;
            else
                allWrong = allWrong + 1;
            end
        end
    end
    
    fprintf('All right: %f\n', allRight/total);
    fprintf('2 of 3 right: %f\n', kindaRight/total);
    fprintf('2 of 3 wrong: %f\n', kindaWrong/total);
    fprintf('All wrong: %f\n', allWrong/total);
    
end

function brainSV = removeBackgroundSV(superVoxels)
    
%threshold can be adjusted based on testing, but since we're only
%trying to remove background, it should be very low intensity
    threshold = 0.03;
    findByIntensity = false;
    
    if findByIntensity
        minIntensity = min(superVoxels(:,10));
        maxIntensity = max(superVoxels(:,10));
        intensityThreshold = minIntensity + (maxIntensity - minIntensity)*threshold;

        
        brainInd = find(superVoxels(:,10) > intensityThreshold);
        brainSV = superVoxels(brainInd,:);
        backgroundInd = find(superVoxels(:,10) < intensityThreshold);
        backgroundSV = superVoxels(backgroundInd,:);
        
        %creating bounding box for finding ventricals
        boxX = [min(brainSV(:,1)),max(brainSV(:,1))];
        boxY = [min(brainSV(:,2)),max(brainSV(:,2))];
        boxZ = [min(brainSV(:,3)),max(brainSV(:,3))];
        
        ventricals = find((boxX(1) < backgroundSV(:,1)) & ...
                          (boxX(2) > backgroundSV(:,1)) & ...
                          (boxY(1) < backgroundSV(:,2)) & ...
                          (boxY(2) > backgroundSV(:,2)) & ...
                          (boxZ(1) < backgroundSV(:,3)) & ...
                          (boxZ(2) > backgroundSV(:,3)));
        if size(ventricals,1) ~= 0
            numBrSV = size(brainSV,1);
            numVentSV = size(ventricals,1);
            numFeat = size(brainSV,2);
            newbrain = zeros(numBrSV + numVentSV, numFeat);
            newbrain(1:numBrSV,:) = brainSV;
            newbrain((numBrSV+1):end,:) = backgroundSV(ventricals,:);
        end
        fprintf('There are %d ventricals SVs in this brain\n', ...
                size(ventricals,1));
        fprintf('There are %d nonvent black SVs in this brain\n',...
                size(backgroundSV,1) - size(ventricals,1));
    else
        brainInd = find((superVoxels(:,7) >= threshold) | ...
                        (superVoxels(:,8) >= threshold) | ...
                        (superVoxels(:,9) >= threshold));
        brainSV = superVoxels(brainInd,:);
        fprintf('%d of %d SV kept\n',size(brainInd,1), ...
                size(superVoxels,1));
    end
    
end

function [AD_SVMstruct CN_SVMstruct MCI_SVMstruct] = trainSVM(training, groups)
    
    global trainMCI
    
    
    % Display can be changed to 'iter' for better readouts during the
    % process and 'final' for the final display once
    % converged. MaxIter is originally 15000, but it seems to scale
    % with the size of the training data, but the tests have been
    % too few to see for certain
    
    options = statset('Display', 'iter', 'MaxIter', 45000);
    
    fprintf('Training AD SVM...\n');
    AD_SVMstruct = svmtrain(training,groups{1},'options',options,...
                            'kernel_function','mlp');
    if trainMCI
        fprintf('Training MCI SVM...\n');
        MCI_SVMstruct = svmtrain(training,groups{2},'options',options,...
                                 'kernel_function','mlp');
        fprintf('Training CN SVM...\n');
        CN_SVMstruct = svmtrain(training,groups{3},'options',options,...
                                'kernel_function','mlp');
    else
        fprintf('Training CN SVM...\n');
        CN_SVMstruct = svmtrain(training,groups{2},'options',options,...
                                'kernel_function','mlp');
    end
end