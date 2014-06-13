%% This file contains code segments that are currently unused.  The
% initial unusedCode function is a dummy so that MATLAB formatting
% while be maintained

function unusedCode()

    disp('Unused code is commented within source file')
    
    %% From CRF_fastTune.m
    % -------- from CRF_fastTune function, after load_nifti, before
    % makeCrossFold
    
    % Leave One Out
    %testing = [leaveOut];
    %if fold == 1
    %    training = 2:nExamples;
    % elseif fold == nExamples
    %     training = 1:(nExamples-1);
    %else
    %    training = [1:(fold-1) (fold + 1):nExamples];
    % end

    % Random split
    % [training, testing] = splitTrainTest(80,20, nExamples); %Make
    % Training/Testing Vectors
    
    
    
    % ---------- from CRF_fastTune function, after Creating Adj Matrix section
    % Calculate Other Features

    %attepted to make binning feature. Never helped. 
    %[xCor, yCor] = cor_feats(nRows, nCols, nNodes); 
    % [WMMin, GMMin, CFMin, BGMin] = min_bins(X,y,nExamples,training);
    % 
    % % Show min values
    % for j = 1:length(testing)
    %     i = testing(j);
    %     Mins = BGMin{i} + (CFMin{i}*2) + (GMMin{i} * 3) + (WMMin{i} * 4);
    %     Mins = reImage(Mins, ZmaskFlat{i});
    %     lenT = length(testing);
    %     lenT = sqrt(lenT);

    %     subplot(lenT, lenT+1, j);
    %     imagesc(reshape(Mins,nRows,nCols));
    %     colormap gray
    % end
    % suptitle('Min Values');
    % if pauses
    %     fprintf('(paused)\n');
    %     pause
    % end;

    
    %-------- CRF_fastTune function, at end of function
        % %% Train with Pseudo-likelihood
    % w = startW;
    % trainingEx = examples(training(:));
    % funObj = @(w)UGM_CRF_PseudoNLL(w,trainingEx(:).Xnode,trainingEx(:).Xedge,y(training),trainingEx(:).nodeMap,trainingEx(:).edgeMap,trainingEx(:).edgeStruct);
    % options.maxFunEvals = 20;
    % 
    % lambda = ones(size(w));
    % penalizedFunObj = @(w)penalizedL2(w,funObj,lambda);
    % 
    % w = minFunc(penalizedFunObj,w);
    % 
    % fprintf('ICM Decoding with estimated parameters...\n');
    % decode(w,testing,testing,nRows,nCols,origY, 'Psudo-likelihood Decoding',ZmaskFlat);   
    % if pauses
    %     fprintf('(paused)\n');
    %     pause
    % end


    % %% Loopy belief propagation
    % w = startW;
    % whos
    % funObj = @(w)UGM_CRFcell_NLL(w,examples,@UGM_Infer_LBP);
    % %w = minConf_TMP(funObj,w,LB,UB);
    % options.maxFunEvals = 6;
    % 
    % %L2 Normalization
    % lambda = ones(size(w));
    % penalizedFunObj = @(w)penalizedL2(w,funObj,lambda);
    % 
    % w = minFunc(penalizedFunObj, w);
    % decode(w,examples,testing,nRows,nCols,y, 'LBP Decoding');   
    % if pauses
    %     fprintf('(paused)\n');
    %     pause
    % end
    % 
end


