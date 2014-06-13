%Matlab paths that need to be added to run:
%addpath(genpath('acmi/fmri/UGM/UGM'))
%addpath(genpath('acmi/fmri/UGM/CRFcell')

%Training begins at 121. The function decode at 219 is where the decoding
%happens. Currently SGD training and ICM with restarts work best and are
%uncommented. Other methods are commented out. 

function CRFGM2_2D()
clear all
close all
pauses = 0; %turn pausing on/off


%% Initialization
[X, y, nExamples] = load_ibsr('/acmi/chris13/analyze_IB/20Normals_T1_brain/', '/acmi/chris13/20Normals_T1_seg/'); %Load Current IBSR files
[training, testing] = splitTrainTest(90,10, nExamples); %Make Training/Testing Vectors

%% Initial Plots
%MRI's
figure;
lenT = length(testing);
lenT = sqrt(lenT);
for i = 1:length(testing)
    subplot(floor(lenT+1), floor(lenT+1), i);
    imagesc(reshape(X{testing(i)}(:,:),256,256)); %TODO: get rid of this last hardcoded size
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
for i = 1:length(testing)
    subplot(floor(lenT+1), floor(lenT+1), i);
    imagesc(reshape(y{testing(i)}(:,:),256,256));
    colormap gray
end
suptitle('Segmentation Truth');
if pauses
    fprintf('(paused)\n');
    pause
end;

%% Make Average Neighbor intensity a feature
nBors = make_nBors(X, nExamples);

%% Make X,y Into Correct Shape and correct Bias

[nRows,nCols] = size(X{1});
nNodes = nRows*nCols;

[origX, origY, Zmask, X, y] = maskZeros(X, y, nExamples);

nStates = max(y{1}(:)); %assume y{1} has all states
ZmaskFlat = cell(20);
for i=1:nExamples
    nPixels = size(X{i},1);
    
    y{i} = reshape(y{i},[1, 1 nPixels]);
    X{i} = reshape(X{i},1,1,nPixels);
    
    nBors{i} = nBors{i}(Zmask{i});
    nBors{i} = reshape(nBors{i},1,1,nPixels);
    
    ZmaskFlat{i} = reshape(Zmask{i}, 1, 1, nNodes);
end

X = cor_bias(X,nExamples);

%% Make edgeStruts
adj = make_adj(nRows, nCols, nNodes);
examples = cell(nExamples);
for i = 1:nExamples
    maskAdj = adj;
    maskAdj = maskAdj(ZmaskFlat{i},:);
    maskAdj = maskAdj(:,ZmaskFlat{i});
    examples{i}.edgeStruct = UGM_makeEdgeStruct(maskAdj,nStates);
    examples{i}.Y = int32(y{i});
end

%% Calculate Other Features

%[xCor, yCor] = cor_feats(nRows, nCols, nNodes); %TODO mask this
[WMMin, GMMin, CFMin, BGMin] = min_bins(X,y,nExamples,training);

% Show min values
for j = 1:length(testing)
    i = testing(j);
    Mins = BGMin{i} + (CFMin{i}*2) + (GMMin{i} * 3) + (WMMin{i} * 4);
    Mins = reImage(Mins, ZmaskFlat{i});
    lenT = length(testing);
    lenT = sqrt(lenT);
    subplot(floor(lenT+1), floor(lenT+1), j);
    imagesc(reshape(Mins,nRows,nCols));
    colormap gray
end
suptitle('Min Values');
if pauses
    fprintf('(paused)\n');
    pause
end;

%% Make Xnode, Xedge, nodeMap, edgeMap, initialize weights
tied = 1;
for i = 1:nExamples
    %M Xnode
    examples{i}.Xnode = [ones(1,1,size(X{i},3)) X{i} UGM_standardizeCols(nBors{i},tied)]; %GMMin{i} WMMin{i} BGMin{i} CFMin{i}]; %UGM_standardizeCols(nBors{i},tied)]; %add feature matricies here
    
    % Make Xedge
    sharedFeatures = [1 0 1]; %needs to reflect number of features
    examples{i}.Xedge = UGM_makeEdgeFeatures(examples{i}.Xnode,examples{i}.edgeStruct.edgeEnds,sharedFeatures(:));    
    [examples{i}.nodeMap examples{i}.edgeMap w] = UGM_makeCRFmaps(examples{i}.Xnode,examples{i}.Xedge,examples{i}.edgeStruct,0,tied,1,1);
    startW = w;
end

%% Stochastic gradient descent training

stepSize = 1e-4;

for iter = 1:30
    i = training(randsample(length(training),1));
    funObj = @(w)UGM_CRF_NLL(w,examples{i}.Xnode,examples{i}.Xedge,y{i},examples{i}.nodeMap,examples{i}.edgeMap,examples{i}.edgeStruct,@UGM_Infer_LBP);
    [f,g] = funObj(w);
    fprintf('Iter = %d of %d (fsub = %f)\n',iter,30,f);
    w = w - stepSize*g;
end

decode(w,examples,testing,nRows,nCols,origY, 'SGD Decoding', ZmaskFlat);   
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
% trainingEx = examples(training(:));
% funObj = @(w)UGM_CRFcell_NLL(w,examples,@UGM_Infer_LBP);
% %w = minConf_TMP(funObj,w,LB,UB);
% options.maxFunEvals = 6;
% 
% %L2 Normalization
% lambda = ones(size(w));
% penalizedFunObj = @(w)penalizedL2(w,funObj,lambda);
% 
% w = minFunc(penalizedFunObj, w);
%  decode(w,examples,testing,nRows,nCols,origY, 'LBP Decoding', ZmaskFlat);
% if pauses
%     fprintf('(paused)\n');
%     pause
% end


end

%% Stats Functions

function M = mcr(L,R)
%L is the labels, R is the results, returns the misclassification rate
dif = L - R;
K = dif(:,:,:) ~= 0;
M = sum(sum(sum(K)));
M = (M/numel(L))*100;
fprintf('\nMCR: %f', M);
end

function M = vOverlap(L,R)
%L is the labels, R is the results, returns the % volume overlap
V1 = sum(sum(sum(L == 4)));
V2 = sum(sum(sum(R == 4)));
[a,b,c] = size(L);
dif = ismember(L,4)+ismember(R,4);
K = dif(:,:,:) == 2;
M = sum(sum(sum(K)));
M = ((M*2)/(abs(V1+V2)))*100;
fprintf('\n%% Volume Overlap: %f', M);
end

function M = vComp(L,R)
V1 = sum(sum(sum(L == 4)));
V2 = sum(sum(sum(R == 4)));
K = abs(V1-V2);
M = sum(sum(sum(K)));
M = (M/(abs(V1+V2)/2))*100;
fprintf('\n%% Volume Difference: %f', M);
end

%% Other Functions
function [train, test] = splitTrainTest(tr,te,number)
all = randperm(number);
nTrain = int32((tr/(tr+te)) * number);
train = all(1:nTrain);
test = all(nTrain:end);
end

function decode(w,examples,testing,nRows,nCols,y,plotTitle, ZmaskFlat)
figure;
lenT = length(testing);
lenT = sqrt(lenT);
aMCR = 0;
aOver = 0;
aDif = 0;
for i = 1:length(testing)
    j =testing(i);
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,examples{j}.Xnode,examples{j}.Xedge,examples{j}.nodeMap,examples{j}.edgeMap,examples{j}.edgeStruct,1);
    
    %Other decoding methods(some not yet working with cell arrays)
    %yDecode = UGM_Decode_LBP(nodePot,edgePot,examples{j}.edgeStruct);
    %yDecode = UGM_Decode_ICM(nodePot,edgePot,examples{i}.edgeStruct);
    %yDecode = int32(UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_LBP));
    %yDecode2 = UGM_Infer_LBP(nodePot,edgePot,examples{j}.edgeStruct);
    yDecode = UGM_Decode_ICMrestart(nodePot,edgePot,examples{j}.edgeStruct,60); %last value is number of restarts
    
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
    subplot(floor(lenT+1),floor(lenT+1), i);
    imagesc(reshape(yDecode,nRows,nCols));
    colormap gray

    %Evaluate
    aMCR = aMCR + mcr((reshape(y{testing(i)},nRows,nCols)),reshape(yDecode,nRows,nCols));
    aOver = aOver + vOverlap((reshape(y{testing(i)},nRows,nCols)),reshape(yDecode,nRows,nCols));
    aDif = aDif + vComp((reshape(y{testing(i)},nRows,nCols)),reshape(yDecode,nRows,nCols));
    
end
 fprintf('\n Average MCR: %f\n Average Volume Overlap: %f\n Average Volume Difference: %f\n', aMCR/length(testing), aOver/length(testing), aDif/length(testing));
suptitle(plotTitle);
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

%% File Loading

function [X,y,nExamples] = load_ibsr(imDir, segDir)
offsets = [0 1 2 1 1 0 0 0 0 3 3 -4 3 2 6 8 1 0 0 2];
%Load data from ibsr database
bList = dir(strcat(imDir,'*.hdr'));
rawImages = cell(length(bList));
nExamples = length(bList);
for f = 1:length(bList)
    fid = fopen(strcat(imDir,bList(f).name));
    stats = fscanf(fid, '%d');
    nSlices = stats(3)+2;
    fid = fopen(strcat(imDir,bList(f).name(1:end-4),'.buchar'));
    I = fread(fid,inf,'uint8=>int32');
    I = reshape(I, 256, 256, nSlices);
    rawImages{f} = reshape(I(:,:,int32(30+offsets(f))),256,256); %just grabbing middle for now
end

bList = dir(strcat(segDir,'*.hdr'));
segs = cell(length(bList));
for f = 1:length(bList)
    fid = fopen(strcat(segDir,bList(f).name));
    stats = fscanf(fid, '%d');
    nSlices = stats(3);
    fid = fopen(strcat(segDir,bList(f).name(1:end-4),'.buchar'));
    I = fread(fid,inf,'uint8=>int32');
    I = reshape(I, 256, 256, nSlices);
    segs{f} = reshape(I(:,:,int32(30)),256,256);
end
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
end

function [X,y,nExamples] = load_nifti()
%not functional yet. Shows how to uncompress and load file, however
%for NIFTI files(needs spm8 in path)
data = ['/acmi/fmri/AD_T1/patient1.nii';'/acmi/fmri/AD_T1/patient3.nii';'/acmi/fmri/AD_T1/patient2.nii'];
flist_names_T1 = cellstr(data);
flist_names_T1{1,1}
namesCount = 0;
I_t1uncompress = wfu_uncompress_nifti('/acmi/fmri/AD_T1/patient1.nii');
I_uncompt1 = spm_vol(I_t1uncompress);
I_T1 = spm_read_vols(I_uncompt1);
I_T1 = (I_T1 - min(I_T1(:)))*10/(max(I_T1(:)) - min(I_T1(:)));

I = (I - min(I(:)))*255/(max(I(:)) - min(I(:))); %Threshold to 255
[a,b,c] = size(I);
X = reshape(I_T1(100,150:160,100:110),11,11);

X = reshape(I(:,:,11),a,b);
y = int32(reshape(y(:,:,8),a,b));
nExamples = 1;
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

function adj = make_adj(nRows, nCols, nNodes)
%Need to test diagonal edges more
adj = sparse(nNodes,nNodes);

% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

% Add Down-Right Edges
ind = 1:nNodes;
exclude = union(sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols), ...
    sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows]))); % No Down edge for last row or Column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1+nRows)) = 1;

% Add Down-Left Edges
ind = 1:nNodes;
exclude = union(sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols), ...
    sub2ind([nRows nCols],1:nRows,repmat(1,[1 nRows]))); % No Down edge for last row or Column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1-nRows)) = 1;

% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;

% Add Up/Left Edges
adj = adj+adj';
end
