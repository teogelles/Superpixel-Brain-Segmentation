%% CRF_test
% Primary author: Chris Magnano
%
% Editors: Teo Gelles
%          Andrew Gilchrist-Scott
%
% This file contains the main script for training the CRF using
% every image in IBSR

%Matlab paths that need to be added to run:
%addpath(genpath('/sonigroup/summer2014/USER/brainseg2014/UGM'))
%addpath(genpath('/sonigroup/fmri/spm8'))


function CRF_train(iterations,res,usePriors,useCdist)
    
    if (usePriors == true && useCdist == true)
        exception = MException('Parameters:Conflicting',['Cannot use both Priors ' ...
                            'and Cdist']);
        throw(exception);
    end
    
    close all
    
    global saveFile tempDir;
    saveFile = makeSaveFile(res, iterations, usePriors, useCdist);
    tempDir = '/scratch/tgelles1/summer2014/training/tempDir/';

    % Load IBSR Data
    [X, y, nExamples] = load_nifti('/sonigroup/fmri/IBSR_nifti_stripped/', res);
    
    training = 1:nExamples;
    
    % Init Features
    % Make Average Neighbor intensity a feature
    fprintf('Creating Neighborhood feature...\n');
    nBors = make_nBors(X, nExamples);
    cDist = NaN;
    if useCdist % If results of spm8 are being used, there is
                % not enough memory on the swatcs computers
                % to use the cDist feature
        
        % Make distance to center a feature
        fprintf('Creating Distance to Center feature...\n')
        cDist = center_distance(X, nExamples);    
    end
    
    
    % Make X,y Into Correct Shape and correct Bias   
    [origX, origY, Zmask, ZmaskFlat, X, y, nStates, sizes, nPixelsArray, ...
     nBors, cDist] = reshapeMatrices(nExamples, X, y, nBors, cDist, ...
                                     useCdist);

    
    % Make edgeStructs COULD PARALLELIZE
    examples = makeEdgeStructs(nExamples, nStates, origX, y, ...
                               sizes, ZmaskFlat);

    clear('y');
    clear('origX');
    clear('sizes')
    
    [examples, w] = prepareExamples(nExamples, examples, res, Zmask, ...
                               nPixelsArray, X, nBors, cDist, ...
                               usePriors, useCdist);
    
    clear('X');
    clear('nBors');
    clear('cDist');
    clear('nPixelsArray');
    clear('Zmask');
    
    % Stochastic gradient descent training
    trainCRF(nExamples, examples, iterations, training, w, origY, ZmaskFlat)
end

function [origX, origY, Zmask, X, y] = maskZeros(X, y, nExamples)
    Zmask = cell(nExamples,1);
    newX = cell(nExamples,1);
    newY = cell(nExamples,1);
    parfor i = 1:nExamples
        Zmask{i} = X{i} ~= 0;
        newX{i} = X{i}(Zmask{i});
        newY{i} = y{i}(Zmask{i});
    end
    origX = X;
    origY = y;
    X = newX;
    y = newY;
end

%% File I/O

function fName = save_data(dir, data, i)
    
    fName = strcat(dir, 'ex', int2str(i));
    
    if ~exist(dir, 'dir')
        mkdir(dir);
    end
    
    save(fName,'data','-v7.3');
end

function loadData = load_data(fName)
    
    loadData = load(fName, 'data');
    loadData = loadData.('data');
end

function [X,y,nExamples] = load_nifti(imDir,res)
%Loads IBSR V2 nifti files

    fprintf('Loading Nifti Images');
    bList = dir(strcat(imDir));
    nExamples = 18; %The IBSR_nifti_stripped directory has 18 images
    X = cell(nExamples,1);
    y = cell(nExamples,1);
    parfor i = 1:nExamples
        
        fprintf('.');
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
    
    fprintf('\n');
    
    parfor i = 1:nExamples
        X{i} = X{i}(1:res:end,1:res:end,1:res:end);
        y{i} = y{i}(1:res:end,1:res:end,1:res:end);
    end
end

%% Feature Making
function nBors = make_nBors(X, nExamples)
%Make Neighbor Intensities another feature
    nBors = cell(nExamples,1);
    rNbor = cell(nExamples,1);
    lNbor = cell(nExamples,1);
    uNbor = cell(nExamples,1);
    dNbor = cell(nExamples,1);

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

function cDist = center_distance(X,nExamples)
    cDist = cell(nExamples,1);
    for i=1:nExamples
        xlen = size(X{i},1);
        ylen = size(X{i},2);
        zlen = size(X{i},3);
        [A,B,C] = meshgrid(1:xlen,1:ylen,1:zlen);
        cDist{i}=sqrt(((A-xlen/2).^2)+((B-ylen/2).^2)+((C-zlen).^2));
    end
    cDist = cor_bias(cDist,nExamples);
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

function adj = make_adj(nRows, nCols, nSlices, nNodes)
    [ii jj] = sparse_adj_matrix([nRows nCols nSlices], 1, inf);
    adj = sparse(ii, jj, ones(1,numel(ii)), nNodes, nNodes);
end

function [ii jj] = sparse_adj_matrix(sz, r, p)
%
% Construct sparse adjacency matrix (provides ii and jj indices into the
% matrix)
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

function priors = load_spm8_matrix(res, tissueNum, imageNum, Zmask, ...
                                   nPixelsArray)
    
    priorsx = load_spm8_priors(res, tissueNum,imageNum);
    
    
    nPixels = nPixelsArray(imageNum);
    
        
    priorsx = priorsx(Zmask{imageNum});
    priorsx = reshape(priorsx,1,1,nPixels);
    
    priors = priorsx;
end

function X = load_spm8_priors(res, tissueNum, imageNum)
% We will only use one of either X or y, but we are unsure of which
% to use as of yet    

    file = '/sonigroup/fmri/IBSR_nifti_stripped/new_Segment_MRF2_dist2/';
    % directory we use for the
    % spm8 tissue segmpentation
    % images

        
    filename = file;
    if (imageNum<10)
        imageNum = strcat('0', int2str(imageNum));
    else
        imageNum = int2str(imageNum);
    end
    
    if (tissueNum == 1)
        filename = strcat(filename, 'rc1IBSR_');
        filename = strcat(filename, imageNum);
        filename = strcat(filename, '_ana.nii');
    elseif (tissueNum == 2)
        filename = strcat(filename, 'rc2IBSR_');
        filename = strcat(filename, imageNum);
        filename = strcat(filename, '_ana.nii');
    elseif (tissueNum == 3)
        filename = strcat(filename, 'rc3IBSR_');
        filename = strcat(filename, imageNum);
        filename = strcat(filename, '_ana.nii');
    end
    
    %Image
    I_t1uncompress = wfu_uncompress_nifti(filename);
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X = I_T1;
    

    %fprintf('%d',size(X,3),size(y,3))
    X = X(1:res:end,1:res:end,1:res:end);
    
end

function saveFile = makeSaveFile(res, numIters, usePriors, useCdist)
    
    global saveFile
    
    saveFile = strcat('/scratch/tgelles1/summer2014/training/', ...
                 'finalWeights-', 'res', int2str(res), 'numIters', ...
                 int2str(numIters));
    
    if (usePriors)
        saveFile = strcat(saveFile, '+spm8priors');
    end
    if (useCdist)
        saveFile = strcat(saveFile, '+centerDistance');
    end
    
    if exist(saveFile, 'file')
        delete(saveFile);
    end
end

function [origX, origY, Zmask, ZmaskFlat, X, y, nStates, sizes, ...
          nPixelsArray, nBors, cDist] = reshapeMatrices(nExamples, ...
                                                      X, y, nBors, ...
                                                      cDist, useCdist);
    global saveFile;
    
    sizes = zeros(nExamples);
    parfor i=1:nExamples
        [nRows,nCols,nSlices] = size(X{i});
        sizes(i) = nRows*nCols*nSlices;
    end

    fprintf('Masking Zeros...\n');
    [origX, origY, Zmask, X, y] = maskZeros(X, y, nExamples);

    nStates = max(y{1}(:)); %assume y{1} has all states 
    ZmaskFlat = cell(nExamples,1);
    
    nPixelsArray = zeros(1, nExamples);

    fprintf('Reshaping Matricies');
    for i=1:nExamples
        fprintf('.');
        nPixels = size(X{i},1);
        nPixelsArray(i) = nPixels;
        
        y{i} = reshape(y{i},[1, 1 nPixels]);
        X{i} = reshape(X{i},1,1,nPixels);
        

        nBors{i} = nBors{i}(Zmask{i});
        nBors{i} = reshape(nBors{i},1,1,nPixels);
        
        if (useCdist)
            cDist{i} = cDist{i}(Zmask{i});
            cDist{i} = reshape(cDist{i},1,1,nPixels);
        end
        
        ZmaskFlat{i} = reshape(Zmask{i}, 1, 1, sizes(i));
    end
    fprintf('\n');
    
    %clear Zmask;
    %we choose to use Zmask later, so it is not cleared here
    fprintf('Correcting Bias...\n');
    X = cor_bias(X,nExamples);
end

function examples = makeEdgeStructs(nExamples, nStates, origX, y, ...
                                    sizes, ZmaskFlat);
    
    global saveFile tempDir;
    
    fprintf('Creating Adj Matrix');
    
    examples = cell(nExamples,1);
    
    for i=1:nExamples
        fprintf('.');
        adj = make_adj(size(origX{i},1), size(origX{i},2), size(origX{i},3),sizes(i));
        adj = adj + adj';
        adj(adj==2) = 1;
        maskAdj = adj;
        maskAdj = maskAdj(ZmaskFlat{i},:);
        maskAdj = maskAdj(:,ZmaskFlat{i});
        clear adj;  %clear things once not needed    
        examples{i}.edgeStruct = UGM_makeEdgeStruct(maskAdj,nStates,1,100);
        clear maskAdj;
        examples{i}.Y = int32(y{i});
        examples{i} = save_data(tempDir, examples{i}, i);
    end

    fprintf('\n');
end

function [examples, w] = prepareExamples(nExamples, examples, res, ...
                                         Zmask, nPixelsArray, X, ...
                                         nBors, cDist, usePriors,useCdist);
    
    global saveFile tempDir;
    
    tied = 1;

    
    fprintf('Creating Xnode, Xedge, and maps');
    
    %Here all the features are put into the final structure
    for i = 1:nExamples
        fprintf('.');
        examples{i} = load_data(examples{i});
        
        %Load spm8 Priors
        if (usePriors)
            % Note: to test space requirements, this code has
            % been changed such that only one priors matrix is loaded
            c1priors = load_spm8_matrix(res, 1, i, Zmask, ...
                                        nPixelsArray);
            %c2priors = load_spm8_matrix(res, 2, i, Zmask, ...
            %                            nPixelsArray);
            %c3priors = load_spm8_matrix(res, 3, i, Zmask, ...
            %                            nPixelsArray);
            
            
            %Make Xnode
            examples{i}.Xnode = [ones(1,1,size(X{i},3)) X{i} ...
                                UGM_standardizeCols(nBors{i},tied), ...
                                c1priors]; %c2priors,
                                           %c3priors];
                                           %GMMin{i}
                                           %WMMin{i}
                                           %BGMin{i}
                                           %CFMin{i}];
                                           %add feature matricies
            sharedFeatures = [1 0 1 0]; % 0 0]; %needs to reflect number
                                        %of features
        elseif (useCdist)
            % We know that useCdist and usePriors are
            % mututally exclusive
            examples{i}.Xnode = [ones(1,1,size(X{i},3)) X{i} ...
                                UGM_standardizeCols(nBors{i},tied) ...
                                cDist{i}];

            %add feature matricies
            sharedFeatures = [1 0 1 1];
            
        else
            examples{i}.Xnode = [ones(1,1,size(X{i},3)) X{i} ...
                                UGM_standardizeCols(nBors{i},tied)];

            %add feature matricies
            sharedFeatures = [1 0 1];
            
        end
        
        examples{i}.Xedge = ...
            UGM_makeEdgeFeatures(examples{i}.Xnode, ...
                                 examples{i}.edgeStruct ...
                                 .edgeEnds,sharedFeatures(:));
        
        %Makes mapping of features to parameters
        [examples{i}.nodeMap examples{i}.edgeMap w] = ...
            UGM_makeCRFmaps(examples{i}.Xnode, ...
                            examples{i}.Xedge, ...
                            examples{i}.edgeStruct,0,tied,1,1);
        
        examples{i} = save_data(tempDir, examples{i}, i);
        
        if (usePriors)
            clear('c1priors')
            %clear('c2priors')
            %clear('c3priors')
        end
    end

    fprintf('\n');
end

function trainCRF(nExamples, examples, iterations, ...
                  training, w, origY, ZmaskFlat);
    n
    global saveFile tempDir;
        
    fprintf('\nBeginning Training\n');
    stepSize = 1e-3;

    %Actual training 
    for iter = 1:iterations
        i = training(randi(length(training),1));
        
        examples{i} = load_data(examples{i});
        %Uncomment to not use mex
        %examples{i}.edgeStruct.useMex = 0;
        %calculate training step
        
        funObj = @(w)UGM_CRF_NLL(w,...
                                 examples{i}.Xnode, ...
                                 examples{i}.Xedge, ...
                                 examples{i}.Y+int32(examples{i}.Y==1), ...
                                 examples{i}.nodeMap, ...
                                 examples{i}.edgeMap, ...
                                 examples{i}.edgeStruct,...
                                 @UGM_Infer_LBP);


        examples{i} = save_data(tempDir, examples{i}, i);
        [f,g] = funObj(w); %calculate gradient from training step
        fprintf('Iter = %d of %d (fsub = %f) on %d\n',iter,iterations,f,i);
        w = w - stepSize*g; %take small step in direction of gradient
    end
    
    save(saveFile, 'w','-v7.3');
    rmdir(tempDir, 's');
end