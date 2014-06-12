%Matlab paths that need to be added to run:
%addpath(genpath('/acmi/chris13/UGM'))
%addpath(genpath('/acmi/fmri/UGM/CRFcell'))
%addpath(genpath('/acmi/fmri/spm8'))


function avg = CRFGM_tune(fold,iterations,leaveOut,res)

close all
pauses = 0; %turn pausing on/off
plots = 0;

fold = fold + 1;

f = 1;
found = 0;
while ~found
    dir = strcat('/local/cmagnan1/tmp',int2str(f),'/'); %directory for temp files
    paramDir = strcat('/local/cmagnan1/params',int2str(f),'/paramsrun'); %directory for parameter files
    if exist(paramDir, 'dir') || exist(dir, 'dir') 
        f = f + 1;
    else
        break;
    end
end
    mkdir(paramDir);

%% TODO's 
% Things to do:
% figure out command line args 
% NEW: make skull distance feature
% MOAR Features

%% Initialization
%[X, y, nExamples] = load_ibsr('/acmi/chris13/analyze_IB/20Normals_T1_brain/', '/acmi/chris13/20Normals_T1_seg/'); %Load IBSR_V1
[X, y, nExamples] = load_nifti('/acmi/fmri/IBSR_nifti_stripped/',res); %Load IBSR_V2
%[X, y, nExamples] = load_ADNI('/acmi/fmri/CN_T1/'); %Load ADNI cog normal


%% Leave One Out
testing = [fold];
if fold == 1
    training = 2:nExamples;
elseif fold == nExamples
    training = 1:(nExamples-1);
else
    training = [1:(fold-1) (fold + 1):nExamples];
end

%% Random split
%[training, testing] = splitTrainTest(80,20, nExamples); %Make Training/Testing Vectors

%% Cross Fold
% nTesting = floor(nExamples * 0.2);
% testStart = 1+nTesting*(fold-1);
% testEnd = (nTesting*fold)+1;
% testing = testStart:testEnd;
% 
% if fold ~=1 || fold ~= 5
%     training = [1:(testStart-1) (testEnd+1):nExamples];
% elseif fold == 1
%     training = (testEnd+1):nExamples;
% else
%     training = 1:testStart-1;
% end


%% Initial Plots
%MRI's
if plots == 1
    figure;
    lenT = length(testing);
    lenT = sqrt(lenT);
    for j = 1:numel(testing)
        subplot(lenT, lenT+1, j);
        size(X{testing(j)})
        i = testing(j);
        imagesc(reshape((X{testing(j)}(:,:,floor(size(X{i},3)/2))),size(X{i},1),size(X{i},2))); %Should not hardcode these
        colormap gray
    end
    suptitle('MRI Images');
    if pauses
        fprintf('(paused)\n');
        pause
    end;
    title('Original MRI');
    
    %Segmentations
    figure;
    lenT = length(testing);
    lenT = sqrt(lenT);
    for j = 1:numel(testing)
        subplot(lenT, lenT+1, j);
        i = testing(j);
        imagesc(reshape((y{testing(j)}(:,:,floor(size(y{i},3)/2))),size(y{i},1),size(y{i},2)));
        colormap gray
    end
    suptitle('Segmentation Truth');
    if pauses
        fprintf('(paused)\n');
        pause
    end
end
%% Make Average Neighbor intensity a feature
fprintf('\nCreating Neighborhood feature...');
nBors = make_nBors(X, nExamples);

%% Make X,y Into Correct Shape and correct Bias

sizes = zeros(nExamples);
for i=1:nExamples
    [nRows,nCols,nSlices] = size(X{i});
    sizes(i) = nRows*nCols*nSlices;
end

fprintf('\nMasking Zeros...');
[origX, origY, Zmask, X, y] = maskZeros(X, y, nExamples);

nStates = max(y{1}(:)); %assume y{1} has all states
ZmaskFlat = cell(nExamples);

fprintf('\nReshaping Matricies');
for i=1:nExamples
    fprintf('.');
    nPixels = size(X{i},1);
    
    y{i} = reshape(y{i},[1, 1 nPixels]);
    X{i} = reshape(X{i},1,1,nPixels);
    
    nBors{i} = nBors{i}(Zmask{i});
    nBors{i} = reshape(nBors{i},1,1,nPixels);
    
    ZmaskFlat{i} = reshape(Zmask{i}, 1, 1, sizes(i));
end

clear Zmask;
fprintf('\nCorrecting Bias...');
X = cor_bias(X,nExamples);

%% Make edgeStruts COULD PARALLELIZE
fprintf('\nCreating Adj Matrix');
examples = cell(nExamples);

for i=1:nExamples
    fprintf('.');
    adj = make_adj(size(origX{i},1), size(origX{i},2), size(origX{i},3),sizes(i)); 
    adj = adj + adj';
    adj(adj==2) = 1;
    maskAdj = adj;
    maskAdj = maskAdj(ZmaskFlat{i},:);
    maskAdj = maskAdj(:,ZmaskFlat{i});
    clear adj;
   % fprintf('\nmaking edge struct on size %d %d %d %d %d\n',size(X{i},1), size(X{i},2), size(X{i},3), sizes(i), numel(X{i}));
    examples{i}.edgeStruct = UGM_makeEdgeStruct(maskAdj,nStates,1,100);%size(X{i},1),size(X{i},2),size(X{i},3));
    clear maskAdj;
    examples{i}.Y = int32(y{i});
    examples{i} = save_data(dir, examples{i}, i);
end

clear y;
%% Calculate Other Features

%[xCor, yCor] = cor_feats(nRows, nCols, nNodes); %TODO mask this
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

%% Make Xnode, Xedge, nodeMap, edgeMap, initialize weights
tied = 1;

origY = save_data(dir, origY, nExamples+1); %save data for later
clear origX;


fprintf('\nCreating Xnode, Xedge, and maps');

for i = 1:nExamples
    fprintf('.');
    examples{i} = load_data(examples{i});
    %Make Xnode
    examples{i}.Xnode = [ones(1,1,size(X{i},3)) X{i} UGM_standardizeCols(nBors{i},tied)]; %GMMin{i} WMMin{i} BGMin{i} CFMin{i}]; %add feature matricies here
    
    % Make Xedge
    sharedFeatures = [1 0 1]; %needs to reflect number of features
    examples{i}.Xedge = UGM_makeEdgeFeatures(examples{i}.Xnode,examples{i}.edgeStruct.edgeEnds,sharedFeatures(:));    
    [examples{i}.nodeMap examples{i}.edgeMap w] = UGM_makeCRFmaps(examples{i}.Xnode,examples{i}.Xedge,examples{i}.edgeStruct,0,tied,1,1);

    examples{i} = save_data(dir, examples{i}, i);
end


%% Stochastic gradient descent training
fprintf('\nBeginning Training\n');
stepSize = 1e-2;

for iter = 1:iterations
    i = leaveOut; %silly way to avoid testing on training
    while i == leaveOut
        i = training(randi(length(training),1));
    end
    examples{i} = load_data(examples{i});
    funObj = @(w)UGM_CRF_NLL(w,examples{i}.Xnode,examples{i}.Xedge,examples{i}.Y+int32(examples{i}.Y==1),examples{i}.nodeMap,examples{i}.edgeMap,examples{i}.edgeStruct,@UGM_Infer_LBP); %this step builds but may benefit from 4 at once
    examples{i} = save_data(dir, examples{i}, i);
    [f,g] = funObj(w);
    fprintf('Iter = %d of %d (fsub = %f) on %d\n',iter,iterations,f,i);
    w = w - stepSize*g; % w is all that is carried between iterations 
    save(paramDir,'w','-v7.3');
end

origY = load_data(origY);

avg = decode(w,examples,testing,origY, 'SGM Decoding', ZmaskFlat, plots, dir);   
if pauses
    fprintf('(paused)\n');
    pause
end


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

%% Stats Functions

function M = mcr(L,R)
%L is the labels, R is the results, returns the misclassification rate
dif = int32(L)-int32(R);
K = dif(:,:,:) ~= 0;
M = sum(sum(sum(K)));
M = (M/numel(R(R~=1)))*100;
fprintf('\nMCR: %f', M);
end

function [results tani] = vOverlap(L,R)
%L is the labels, R is the results, returns the % volume overlap
results = [0 0 0];
tani = [0 0 0];
for i=2:max(max(max(R)))
    max(max(max(R==i)))
    max(max(max(L==i)))
    V1 = sum(sum(sum(L == i)));
    V2 = sum(sum(sum(R == i)));
    dif = ismember(L,i)+ismember(R,i);
    K = dif(:,:,:) == 2;
    N = sum(sum(sum(K)));
    M =((N*2)/(abs(V1+V2)))*100;
    results(i-1) = M;
    tani(i-1) =(N/(V1+V2-N));
    fprintf('\n%% Volume Overlap: %f, Tanimoto: %f', M, N/(V1+V2-N));
end
end

function M = vComp(L,R)
N = 0;
V1 = 0;
V2 = 0;
for i=2:max(max(max(R)))
    V1 = sum(sum(sum(L == i)));
    V2 = sum(sum(sum(R == i)));
    K = abs(V1-V2);
    N = sum(sum(sum(K)));
    M = (N/(abs(V1+V2)/2))*100;
    fprintf('\n%% Volume Difference: %f', M);
end
M = (N/(abs(V1+V2)/2))*100;

end

%% Other Functions
function [train, test] = splitTrainTest(tr,te,number)
all = randperm(number);
nTrain = int32((tr/(tr+te)) * number);
train = all(1:nTrain);
test = all(nTrain:end);
end

function avg = decode(w,examples,testing,y,plotTitle, ZmaskFlat, plots, dir)
if plots == 1
    figure;
end
lenT = length(testing);
lenT = sqrt(lenT);
aGM = 0;
aWM = 0;
aCSF = 0;
aMCR = 0;
aDif = 0;
tWM = 0;
tGM= 0;
tCSF = 0;

for i = 1:length(testing)

    j =testing(i);
    examples{j} = load_data(examples{j});
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,examples{j}.Xnode,examples{j}.Xedge,examples{j}.nodeMap,examples{j}.edgeMap,examples{j}.edgeStruct,1);
    
    %yDecode = UGM_Decode_LBP(nodePot,edgePot,examples{j}.edgeStruct);
    %yDecode = UGM_Decode_ICM(nodePot,edgePot,examples{i}.edgeStruct);
    %yDecode = int32(UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_LBP));
    %yDecode2 = UGM_Infer_LBP(nodePot,edgePot,examples{j}.edgeStruct);
    yDecode = UGM_Decode_ICMrestart(nodePot,edgePot,examples{j}.edgeStruct,30); %last value is number of restarts
    size(yDecode(yDecode==3))
    yDecode = reImage(yDecode, ZmaskFlat{j});
    yDecode(yDecode == 0) = 1;
    
    %plot marginals
    %     yDecode2 = UGM_Infer_LBP(nodePot,edgePot,examples{j}.edgeStruct);
    %     figure;
    %     for k = 1:4
    %         subplot(2,3,k);
    %         marg = reImage(yDecode2(:,k),ZmaskFlat{j});
    %         imagesc(reshape(marg,nRows,nCols));
    %         colormap gray
    %     end
    
    [nRows, nCols, nSlices] = size(y{testing(i)});
    yDecode = reshape(yDecode, nRows, nCols, nSlices);

    if plots == 1
        subplot(lenT,lenT+1, i);
        imagesc(reshape(yDecode(:,:,3),nRows,nCols));
        colormap gray
    end
    size(yDecode)
    
    %Evaluate
    
    aMCR = aMCR + mcr((reshape(y{testing(i)},nRows,nCols,nSlices)),reshape(yDecode,nRows,nCols,nSlices));
    results = [0 0 0];
    [results, tani] =  vOverlap((reshape(y{testing(i)},nRows,nCols,nSlices)),reshape(yDecode,nRows,nCols,nSlices));
    aCSF = aCSF + results(1);
    aGM = aGM + results(2);
    aWM = aWM + results(3);
    tWM = tWM + tani(3);
    tGM = tGM + tani(2);
    tCSF = tCSF + tani(1);
    aDif = aDif + vComp((reshape(y{testing(i)},nRows,nCols,nSlices)),reshape(yDecode,nRows,nCols,nSlices));
    examples{j} = save_data(dir, examples{j}, j);
    
end
 fprintf('\n Average MCR: %f\n Average Volume Overlap: %f %f %f\n Average Volume Difference: %f\nAverage Tanimoto %f %f %f\n', aMCR/length(testing), aCSF/length(testing), ...
    aGM/length(testing), aWM/length(testing), aDif/length(testing), tWM/length(testing), tGM/length(testing), tCSF/length(testing));
 cleanup(dir);
 avg = (aWM + aGM)/(2*length(testing));
%suptitle(plotTitle);
end

function [origX, origY, Zmask, X, y] = maskZeros(X, y, nExamples)
Zmask = cell(nExamples);
newX = cell(nExamples);
newY = cell(nExamples);
for i = 1:nExamples
    Zmask{i} = X{i} ~= 0;
    newX{i} = X{i}(Zmask{i});
    newY{i} = y{i}(Zmask{i});
end
origX = X;
origY = y;
X = newX;
y = newY;
end

function nImage = reImage(masked, ZmaskFlat)
nImage = int32(ZmaskFlat);
nImage(:,:,ZmaskFlat) = masked(:,:,:);
end

function M = detect_skull(I)

[h,w,r] = size(I);

for iCount = 1:r,
    J = imfill(I(:,:,iCount),'holes');
    K = im2bw(J/max(J(:)), 0.3*graythresh(J(:)/max(J(:))));
    [L1,N] = bwlabel(K);
    maxa = 0; maxj=0;
    for jCount=1:N,
        a = sum(sum(L1==jCount));
        if a>maxa,
            maxa=a;
            maxj=jCount;
        end
    end
    L(:,:,iCount) = double(L1==maxj);
end
M = L;
end

%% File I/O

function fName = save_data(dir, data, i)
fName = strcat(dir, 'ex', int2str(i));
if ~exist(dir, 'dir')
    mkdir(dir);
end
save(fName,'data','-v7.3');
% f = fopen(fName, 'w');
% fwrite(f, data, saveClass);
% fclose(f);
end

function loadData = load_data(fName)
% f = fopen(fName, 'r');
% data = fread(f, inf, loadClass);
% fclose(f);
loadData = load(fName, 'data');
loadData = loadData.('data');
end

function cleanup(dir)
if exist(dir, 'dir')
    rmdir(dir, 's');
end
end

function [X,y,nExamples] = load_ibsr_Ana(imDir, segDir)
%% load data for ibsr data type
%needs to take 
bList = dir(strcat(imDir,'*.hdr'));
rawImages = cell(length(bList));
nExamples = length(bList);
%rawImages = zeros(length(bList), 256,256);
for f = 1:length(bList)
    fid = fopen(strcat(imDir,bList(f).name));
    stats = fscanf(fid, '%d');
    %nSlices = stats(3)+2;
    fid = fopen(strcat(imDir,bList(f).name(1:end-4),'.img'));
    I = fread(fid,inf,'uint16=>int32');
    nSlices = size(I,1)/(256*256);
    I = reshape(I, 256, 256, nSlices);
    %rawImages(f,:,:) = reshape(I(:,:,int32(nSlices/2)),256,256);
    rawImages{f} = reshape(I(:,:,int32(nSlices/3)),256,256);
end

bList = dir(strcat(segDir,'*seg_ana.hdr'));
segs = cell(length(bList));
%segs = zeros(length(bList), 256,256);
for f = 1:length(bList)
    fid = fopen(strcat(segDir,bList(f).name));
    %stats = fscanf(fid, '%d');
    %nSlices = stats(3);
    fid = fopen(strcat(segDir,bList(f).name(1:end-4),'.img'));
    I = fread(fid,inf,'uint16=>int32');
    nSlices = size(I,1)/(256*256);
    I = reshape(I, 256, 256, nSlices);
    %segs(f,:,:) = reshape(I(:,:,int32(nSlices/2)),256,256); %Just grab middle for now
    segs{f} = reshape(I(:,:,int32(nSlices/3)),256,256);
end
fclose('all'); %python needs this

%Create X and Y
y = segs;
X = rawImages;
for i=1:nExamples
    %1 is background, 2 is CSF, 3 is GM, 4 is WM
    y{i} = y{i} + 1;%((y{i}==0) + ((y{i}==128)*2) + ((y{i}==192)*3) + ((y{i}==254)*4)); 
    X{i} = double(X{i});
    y{i} = int32(y{i});
        max(y{i}(:))
end
end

function [X,y,nExamples] = load_ibsr(imDir, segDir)
offsets = [0 1 2 1 1 0 0 0 0 3 3 -4 3 2 6 8 1 0 0 2];
%Load data from ibsr database
bList = dir(strcat(imDir,'*.hdr'));
rawImages = cell(length(bList));
nExamples = length(bList)-1; %CAUSE WHAT THE HELL -4



bList = dir(strcat(segDir,'*.hdr'));
segs = cell(length(bList));
ends = zeros(nExamples,1);

for f = 1:length(bList)
    fid = fopen(strcat(segDir,bList(f).name));
    stats = fscanf(fid, '%d');
    nSlices = stats(3);
    ends(f) = nSlices;
    fid = fopen(strcat(segDir,bList(f).name(1:end-4),'.buchar'));
    I = fread(fid,inf,'uint8=>int32');
    I = reshape(I, 256, 256, nSlices);
    if offsets(f) > -1
        segs{f} = I;
    end
end

for f = 1:length(bList)
    fid = fopen(strcat(imDir,bList(f).name));
    stats = fscanf(fid, '%d');
    nSlices = stats(3)+2;
    fid = fopen(strcat(imDir,bList(f).name(1:end-4),'.buchar'));
    I = fread(fid,inf,'uint8=>int32');
    I = reshape(I, 256, 256, nSlices);
    if (ends(f) + offsets(f)) > nSlices
        fprintf('herr!, %d, %d\n', nSlices, ends(f));
        a = segs{f};
        segs{f} = a(:,:,1:nSlices-offsets(f));
        rawImages{f} = I(:,:,offsets(f)+1:end);
    else
        if offsets(f) > -1
            rawImages{f} = I(:,:,offsets(f)+1:ends(f)+offsets(f));
        end
    end
end

segs{12} = segs{11};
rawImages{12} = rawImages{11};

fclose('all'); %python needs this

%Create X and Y
y = segs;
X = rawImages;
for i=1:nExamples
    %1 is background, 2 is CSF, 3 is GM, 4 is WM
    y{i} = ((y{i}==0) + ((y{i}==128)*2) + ((y{i}==192)*3) + ((y{i}==254)*4)); 
    X{i} = double(X{i});
    y{i} = int32(y{i});
end
for i = 1:nExamples
    fprintf('%d, %d, %d\n',size(X{i},3),size(y{i},3), offsets(i))
    X{i} = X{i}(1:2:end,1:4:end,1:2:end);
     y{i} = y{i}(1:2:end,1:4:end,1:2:end);
end
% nExamples = 7;
% X(8:end) = [];
% y(8:end) = [];
end

function [X,y,nExamples] = load_ADNI(imDir)

%load WM
%load GM
%load CSF

bList = dir(strcat(imDir, 'c3*'));

nExamples = length(bList);
nExamples = 4;
X = cell(nExamples,1);
y = cell(nExamples,1);
fprintf('\nLoading %d files',nExamples);
for i = 1:nExamples
    fprintf('.');
    heads = {'m','c1','c2','c3'};
    results = cell(4);
    
    for j=1:4
        heads{j} = strcat(imDir, heads{j}, 'patient', int2str(i), '.nii');
        I_t1uncompress = wfu_uncompress_nifti(heads{j});
        I_uncompt1 = spm_vol(I_t1uncompress);
        I_T1 = spm_read_vols(I_uncompt1);
        results{j} = int32(I_T1);
    end

X{i} = results{1};
marginals = zeros([3 size(results{4})]);
marginals(1,:,:,:) = results{2}(:,:,:);
marginals(2,:,:,:) = results{3}(:,:,:);
marginals(3,:,:,:) = results{4}(:,:,:);

[C y{i}] = max(marginals, [], 1);
y{i} = reshape(((int32(y{i}).*int32(C > 0.5)) + 1), size(X{i})); %Threshold for marginals 

end
fprintf('\n');
for i = 1:nExamples
   % fprintf('%d, %d\n',size(X{i},3),size(y{i},3))
    X{i} = double(X{i}(1:2:end,1:2:end,1:2:end));
     y{i} = double(y{i}(1:2:end,1:2:end,1:2:end));
end
end

function [X,y,nExamples] = load_nifti(imDir,res)
%Loads IBSR V2 nifti files

bList = dir(strcat(imDir));
nExamples = 18; %for README
X = cell(nExamples);
y = cell(nExamples);
for i = 1:nExamples
    
    if i < 10
        place = strcat(imDir,'IBSR_0',int2str(i),'/');
        fileHead = strcat('IBSR_0',int2str(i));
    else
        place = strcat(imDir,'IBSR_',int2str(i),'/');
        fileHead = strcat('IBSR_',int2str(i));
    end
    
    %Image
    I_t1uncompress = wfu_uncompress_nifti(strcat(place,fileHead,'_ana_strip.nii'));
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X{i} = I_T1;
    
    %Segmentation
    I_t1uncompress = wfu_uncompress_nifti(strcat(place,fileHead,'_segTRI_ana.nii'));
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    y{i} = I_T1+1;
    
    
end
for i = 1:nExamples
    fprintf('%d, %d\n',size(X{i},3),size(y{i},3))
    X{i} = X{i}(1:res:end,1:res:end,1:res:end);
     y{i} = y{i}(1:res:end,1:res:end,1:res:end);
end
end

function [X,y,nExamples] = load_brainweb()
% Currently Only Works for 1 brain. 
if exist('/home/cmagnan1/phantom_full.rawb', 'file')
  fid = fopen('/home/cmagnan1/phantom_full.rawb');
  I = fread(fid,inf,'uint8=>int32');
  I = reshape(I, 181, 217, 181);
end

if exist('/home/cmagnan1/phantom_1.0mm_normal_bck.rawb', 'file')
  fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_bck.rawb');
  bck = fread(fid,inf,'uint8=>int32');
  bck = reshape(bck, 181, 217, 181);
end

if exist('/home/cmagnan1/phantom_1.0mm_normal_csf.rawb', 'file')
  fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_csf.rawb');
  csf = fread(fid,inf,'uint8=>int32');
  csf = reshape(csf, 181, 217, 181);
end

if exist('/home/cmagnan1/phantom_1.0mm_normal_gry.rawb', 'file')
  fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_gry.rawb');
  gry = fread(fid,inf,'uint8=>int32');
  gry = reshape(gry, 181, 217, 181);
end

if exist('/home/cmagnan1/phantom_1.0mm_normal_wht.rawb', 'file')
  fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_wht.rawb');
  wht = fread(fid,inf,'uint8=>int32');
  wht = reshape(wht, 181, 217, 181);
end


%Create label data
bck = (bck>=30) * 1;
csf = (csf>=30) * 2;
wht = (wht>=30) * 3;
gry = (gry>=30) * 4;
y = bck + wht;
y(y==3) = 2;
y = y + gry;
y(y>3) = 3;
y = y + csf;
y(y>4) = 4;

X = I;
nExamples = 1;
end

%% Feature Making
function nBors = make_nBors(X, nExamples)
%Make Neighbor Intensities another feature
nBors = cell(nExamples);
rNbor = cell(nExamples);
lNbor = cell(nExamples);
uNbor = cell(nExamples);
dNbor = cell(nExamples);

for i = 1:nExamples
    rNbor{i} = zeros(size(X{i}));
    lNbor{i} = zeros(size(X{i}));
    uNbor{i} = zeros(size(X{i}));
    dNbor{i} = zeros(size(X{i}));
    
    rNbor{i}(1:end-1,:) = rNbor{i}(1:end-1,:) + X{i}(2:end,:);
    lNbor{i}(2:end,:) = lNbor{i}(2:end,:) + X{i}(1:end-1,:);
    uNbor{i}(:,1:end-1) = uNbor{i}(:,1:end-1) + X{i}(:,2:end);
    dNbor{i}(:,2:end) = dNbor{i}(:,2:end) + X{i}(:,1:end-1);
    
    nBors{i} = rNbor{i} + lNbor{i} + uNbor{i} + dNbor{i};
    nBors{i} = nBors{i} ./ 4;
end
end

function X = cor_bias(X, nExamples)
%Intensity Bias Correction
for i=1:nExamples
%     notZ = X{i}((X{i} ~= 0)); 
%     ignoreStd = std(notZ);
%     ignoreMean = mean(notZ);
%     X{i} = (X{i}-ignoreMean)/ignoreStd + 4;
    X{i} = UGM_standardizeCols(X{i},1);
end
end

function [WMMin, GMMin, CFMin, BGMin] = min_bins(X, y, nExamples, training)
GM = 0;
WM = 0;
CF = 0;
BG = 0;
for i = 1:length(training)
    GMbin = (y{training(i)} == 3);
    GM = GM + sum(sum(X{training(i)} .* GMbin))/sum(sum(GMbin)); %average GM value
    
    WMbin = (y{training(i)} == 4);
    WM = WM + sum(sum(X{training(i)} .* WMbin))/sum(sum(WMbin));
    
    BGbin = (y{training(i)} == 1);
    if sum(sum(BGbin == 0))
        BG = BG + -5;
    else
        BG = BG + sum(sum(X{training(i)}.* BGbin))/sum(sum(BGbin));
    end
    
    CFbin = (y{training(i)} == 2);
    CF = CF + sum(sum(X{training(i)} .* CFbin))/sum(sum(CFbin));
end

WM = WM/length(training);
GM = GM/length(training) - 1;
CF = CF/length(training) - 1;
BG = BG/length(training);

fprintf('Average Values: %f, %f, %f, %f \n', WM, GM, BG, CF);

BGMin = cell(nExamples);
WMMin = cell(nExamples);
GMMin = cell(nExamples);
CFMin = cell(nExamples);
p = 1;
figure;
for i = 1:nExamples
    WMDif = abs(X{i} - WM);
    GMDif = abs(X{i} - GM);
    BGDif = abs(X{i} - BG);
    CFDif = abs(X{i} - CF);
    
    BGMin{i} = BGDif < (((GMDif < WMDif) .* GMDif) + ((WMDif <= GMDif) .* WMDif));
    WMMin{i} = WMDif <= (((GMDif < BGDif) .* GMDif) + ((BGDif <= GMDif) .* BGDif));
    GMMin{i} = GMDif < (((BGDif < WMDif) .* BGDif) + ((WMDif <= BGDif) .* WMDif));
    CFMin{i} = CFDif < (((GMDif < WMDif) .* GMDif) + ((WMDif <= GMDif) .* WMDif));
    
    BGMin{i} = BGMin{i} .* (BGDif <= CFDif);
    WMMin{i} = WMMin{i} .* (WMDif <= CFDif);
    GMMin{i} = GMMin{i} .* (GMDif <= CFDif);
    CFMin{i} = CFMin{i} .* (CFDif <= BGDif);
end
end

function [xCor yCor] = cor_feats(nRows, nCols, nNodes)
%% Make x and y positions another feature
xCor = 1:nCols;
yCor = 1:nRows;
xCor = reshape(xCor, nCols, 1);
yCor = reshape(yCor, nRows, 1);
xCor = reshape(repmat(xCor, [1 nRows]),nRows,nCols);
yCor = repmat(yCor, [1 nCols]);
xCor = double(reshape(xCor,1,1,nNodes));
yCor = double(reshape(yCor,1,1,nNodes)); 
xCor = abs(xCor - nCols/2);
yCor = abs(yCor - nRows/2);

yCor = repmat(yCor, [1 1 1]);
xCor = repmat(xCor, [1 1 1]);
%cDist = sqrt(xCor.^2 + yCor.^2); %change to single distance from center\
end

function adj = make_adj(nRows, nCols, nSlices, nNodes)
[ii jj] = sparse_adj_matrix([nRows nCols nSlices], 1, inf);
adj = sparse(ii, jj, ones(1,numel(ii)), nNodes, nNodes);
% %Need to test diagonal edges more
% adj = sparse(nNodes,nNodes);
% 
% % Add Down Edges
% ind = 1:nNodes-1;
% exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
% ind = setdiff(ind,exclude);
% adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;
% 
% % Add Down-Right Edges
% ind = 1:nNodes-1;
% exclude = union(sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols), ...
%     sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows]))); % No Down edge for last row or Column
% ind = setdiff(ind,exclude);
% adj(sub2ind([nNodes nNodes],ind,ind+1+nRows)) = 1;
% 
% 
% %adj(sub2ind([nNodes nNodes], ind, ind+(nRows*nCols))) = 1;
% 
% % Add Down-Left Edges
% ind = 1:nNodes-1;
% exclude = union(sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols), ...
%     sub2ind([nRows nCols],1:nRows,repmat(1,[1 nRows]))); % No Down edge for last row or Column
% ind = setdiff(ind,exclude);
% adj(sub2ind([nNodes nNodes],ind,ind+1-nRows)) = 1;
% 
% % Add Right Edges
% ind = 1:nNodes-1;
% exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
% ind = setdiff(ind,exclude);
% adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;
% 
% % Add Up/Left Edges
% adj = adj+adj';
end

function [ii jj] = sparse_adj_matrix(sz, r, p)
%
% Construct sparse adjacency matrix (provides ii and jj indices into the
% matrix)
%
% Usage:
%   [ii jj] = sparse_adj_matrix(sz, r, p)
%
% inputs:
%   sz - grid size (determine the number of variables n=prod(sz), and the
%        geometry)
%   r  - the radius around each point for which edges are formed
%   p  - in what p-norm to measure the r-ball, can be 1,2 or 'inf'
%
% outputs
%   ii, jj - linear indices into adjacency matrix (for each pair (m,n)
%   there is also the pair (n,m))
%
% How to construct the adjacency matrix?
% >> A = sparse(ii, jj, ones(1,numel(ii)), prod(sz), prod(sz));
%
%
% Example:
% >> [ii jj] = sparse_adj_matrix([10 20], 1, inf);
% construct indices for 200x200 adjacency matrix for 8-connect graph over a
% grid of 10x20 nodes.
% To visualize the graph:
% >> [r c]=ndgrid(1:10,1:20);
% >> A = sparse(ii, jj, ones(1,numel(ii)), 200, 200);;
% >> gplot(A, [r(:) c(:)]);
%
%
%
% Copyright (c) Bagon Shai
% Department of Computer Science and Applied Mathmatics
% Wiezmann Institute of Science
% http://www.wisdom.weizmann.ac.il/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% Sep. 2010
%

% number of variables
n = prod(sz);
% number of dimensions
ndim = numel(sz);

tovec = @(x) x(:);
N=cell(ndim,1);
I=cell(ndim,1);

% construct the neighborhood
fr=floor(r);
for di=1:ndim
    N{di}=-fr:fr;
    I{di}=1:sz(di);
end

[N{1:ndim}]=ndgrid(N{:});
[I{1:ndim}]=ndgrid(I{:});
N = cellfun(tovec, N, 'UniformOutput',false);
N=[N{:}];
I = cellfun(tovec, I, 'UniformOutput',false);
I=[I{:}];
% compute N radius according to p
switch lower(p)
    case {'1','l1',1}
        R = sum(abs(N),2);
    case {'2','l2',2}
        R = sum(N.*N,2);
        r=r*r;
    case {'inf',inf}
        R = max(abs(N),[],2);
    otherwise
        error('sparse_adj_matrix:norm_type','Unknown norm p (should be either 1,2 or inf');
end
N = N(R<=r+eps,:);

% "to"-index (not linear indices)
ti = bsxfun(@plus, permute(I,[1 3 2]), permute(N, [3 1 2]));
sel = all(ti >= 1, 3) & all( bsxfun(@le, ti, permute(sz, [1 3 2])), 3);
csz = cumprod([1 sz(1:(ndim-1))]);
jj = sum( bsxfun(@times, ti-1, permute(csz, [1 3 2])), 3)+1; % convert to linear indices
ii = repmat( (1:n)', [1 size(jj,2)]);
jj = jj(sel(:));

ii = ii(sel(:));

end
