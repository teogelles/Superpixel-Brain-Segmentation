%% getADNITissues
% This program uses a CRF to determine the tissue segmentation on
% ~300 ADNI brain images.
%
% Written by Teo Gelles and Andrew Gilchrist-Scott
% Advised by Ameet Soni

function getADNITissues(type,section,weightFile,useCdist,usePriors)
% @param type - Type of subject (AD,MCI, or CN)
% @param section - portion of each types data we want to run. This
% is useful for using condor on this functoin sinc we can decode
% multiple images simultaneously. Should be in [1,10]
% @param weightFile - filename of .mat file with the weights from
% the optimal CRF training using the full 18 IBSR images
    
% Setting default parameters
    if ~exist('weightFile','var')
        weightFile = ['/scratch/tgelles1/summer2014/training/' ...
                      'finalWeights-res1numIters280+' ...
                      'centerDistance.mat'];
        useCdist = true;
        usePriors = false;
    end    
    if ~exist('useCdist','var')
        useCdist = true;
    end
    
    if ~exist('usePriors','var')
        usePriors = false;
    end
    
    checkInput(type,section,weightFile);
    
    global dir;
    dir = makeDirs(usePriors, useCdist);
        
    % Load the data
    weights = load(weightFile);
    weights = weights.w;
    [images indImages] = load_nifti(type,section);
    %%% For testing, setting numImages to 1 for fast runs
    numImages = 1; %size(images,1);
    
    fprintf('Creating Neighborhood feature...\n');
    nBors = make_nBors(images, numImages);
    cDist = NaN;
    if useCdist % If results of spm8 are being used, there is
                % not enough memory on the swatcs computers
                % to use the cDist feature
        
        % Make distance to center a feature
        fprintf('Creating Distance to Center feature...\n')
        cDist = center_distance(images, numImages);    
    end
    
    % Make images,y Into Correct Shape and correct Bias   
    [origimages, Zmask, ZmaskFlat, images, nStates, sizes, nPixelsArray, ...
     nBors, cDist] = reshapeMatrices(numImages, images, nBors, cDist, ...
                                     useCdist);
    
    brainCRF = makeEdgeStructs(numImages, nStates, origimages, ...
                               sizes, ZmaskFlat);
    
    clear sizes;

    brainCRF = prepareBrainCRF(numImages, brainCRF, Zmask, ...
                               nPixelsArray, images, nBors, cDist, ...
                               usePriors, useCdist);
    
    decode(weights,brainCRF,ZmaskFlat,origimages,numImages,indImages);
    
    fprintf('Done\n');
end

function checkInput(type,section,weightfile)
% Checks if the inputs are valid
    
    if ~(strcmp(type,'AD') || strcmp(type,'MCI') || strcmp(type, ...
                                                          'CN'))
        excep = MException('ArgumentError:type',['Type %s is not ' ...
                            'valid. Please enter AD, MCI, or CN.'], ...
                           type);
        throw(excep);
    end
    
    if (section < 1) || (section > 10)
        excep = MException('ArgumentError:section',['Section %d is ' ...
                            'inappropriate. Please input a section ' ...
                            'number between 1 and 10.']);
        throw(excep);
    end
    
    if ~exist(weightfile,'file')
        excep = MException('ArgumentError:weightfile',['Weighfile ' ...
                            'does not exist or you do not have the ' ...
                            'permissions to access it.']);
        throw(excep);
    end
end

function [images patients] = load_nifti(type,section)
%Loads IBSR V2 nifti files

    fprintf('Loading Nifti Images');
    
    patients = getPatients(type,section);
    
    images = cell(size(patients,2),1);
    
    
    parfor pat_i = patients
        
        fprintf('.');
        filename = strcat('/acmi/fmri/ADNI_Stripped/',type,...
                          sprintf('%03d',pat_i),'.nii');
        
        %Image
        I_t1uncompress = wfu_uncompress_nifti(filename);
        I_uncompt1 = spm_vol(I_t1uncompress);
        I_T1 = spm_read_vols(I_uncompt1);
        images{pat_i} = I_T1;
    end
    
    fprintf('\n');
end

function patients = getPatients(type, section)
    
    totalPatients = [92, 203 ,102];
    types = {'AD','MCI','CN'};
    for i = 1:3
        if strcmp(type,types{i})
            patientsOfType = totalPatients(i);
            break
        end
    end
    
    step = ceil(patientsOfType/10);
    if section == 1
        patients = 1:step;
    elseif section == 10
        patients = ((step*9)+1):patientsOfType;
    else
        patients = (section*(step-1) + 1):(section*step);
    end
    patients = patients;
end

function nBors = make_nBors(images, numImages)
%Make Neighbor Intensities another feature
    nBors = cell(numImages,1);
    prNbor = cell(numImages,1);
    lNbor = cell(numImages,1);
    uNbor = cell(numImages,1);
    dNbor = cell(numImages,1);

    for i = 1:numImages
        rNbor{i} = zeros(size(images{i}));
        lNbor{i} = zeros(size(images{i}));
        uNbor{i} = zeros(size(images{i}));
        dNbor{i} = zeros(size(images{i}));
        
        rNbor{i}(1:end-1,:) = rNbor{i}(1:end-1,:) + images{i}(2:end,:);
        lNbor{i}(2:end,:) = lNbor{i}(2:end,:) + images{i}(1:end-1,:);
        uNbor{i}(:,1:end-1) = uNbor{i}(:,1:end-1) + images{i}(:,2:end);
        dNbor{i}(:,2:end) = dNbor{i}(:,2:end) + images{i}(:,1:end-1);
        
        nBors{i} = rNbor{i} + lNbor{i} + uNbor{i} + dNbor{i};
        nBors{i} = nBors{i} ./ 4;
    end
end

function [origimages, Zmask, ZmaskFlat, images, nStates, sizes, ...
          nPixelsArray, nBors, cDist] = reshapeMatrices(numImages, images, nBors, cDist, ...
                                     useCdist);
    
    global dir;
    
    sizes = zeros(numImages);
    parfor i=1:numImages
        [nRows,nCols,nSlices] = size(images{i});
        sizes(i) = nRows*nCols*nSlices;
    end

    fprintf('Masking Zeros...\n');
    [origimages, Zmask, images] = maskZeros(images, numImages);

    nStates = 4; % hardcoded max number of states that we know 
    ZmaskFlat = cell(numImages,1);
    
    nPixelsArray = zeros(1, numImages);

    fprintf('Reshaping Matricies');
    for i=1:numImages
        fprintf('.');
        nPixels = size(images{i},1);
        nPixelsArray(i) = nPixels;
                
        images{i} = reshape(images{i},1,1,nPixels);
        
        
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
    images = cor_bias(images,numImages);
end

function dir = makeDirs(usePriors,useCdist)
    dir = strcat('/scratch/tgelles1/summer2014/ADNI_temp/decodeCD', ...
                 num2str(useCdist),'PR',num2str(usePriors),'/');
    if ~exist(dir,'dir')
        mkdir(dir);
    end
end

function [origimages, Zmask, images] = maskZeros(images, numImages)
    Zmask = cell(numImages,1);
    newimages = cell(numImages,1);
    parfor i = 1:numImages
        Zmask{i} = images{i} ~= 0;
        newimages{i} = images{i}(Zmask{i});
    end
    origimages = images;
    images = newimages;
end


function images = cor_bias(images, numImages)
%Intensity Bias Correction
    for i=1:numImages
        images{i} = UGM_standardizeCols(images{i},1);
    end
end

function brainCRF = makeEdgeStructs(numImages, nStates, origimages, ...
                                    sizes, ZmaskFlat);
    
    global dir;
    
    fprintf('Creating Adj Matrix');
    
    brainCRF = cell(numImages,1);
    
    for i=1:numImages
        fprintf('.');
        adj = make_adj(size(origimages{i},1), size(origimages{i},2), size(origimages{i},3),sizes(i));
        adj = adj + adj';
        adj(adj==2) = 1;
        maskAdj = adj;
        maskAdj = maskAdj(ZmaskFlat{i},:);
        maskAdj = maskAdj(:,ZmaskFlat{i});
        clear adj;  %clear things once not needed    
                    % Unclear why iterations for the below is 100
        brainCRF{i}.edgeStruct = UGM_makeEdgeStruct(maskAdj,nStates,1,100);
        clear maskAdj;
        brainCRF{i} = save_data(dir, brainCRF{i}, i);
    end
    fprintf('\n');
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

function cleanup(dir)
    if exist(dir, 'dir')
        rmdir(dir, 's');
    end
end

function brainCRF = prepareBrainCRF(numImages, brainCRF, ...
                                    Zmask, nPixelsArray, images, ...
                                    nBors, cDist, usePriors,useCdist);
    
    global dir;
    
    tied = 1;
    
    fprintf('Creating imagesnode, imagesedge, and maps');
    
    %Here all the features are put into the final structure
    for i = 1:numImages
        fprintf('.');
        brainCRF{i} = load_data(brainCRF{i});
        
        %Load spm8 Priors
        if (usePriors)
            c1priors = load_spm8_matrix(1, 1, i, Zmask, ...
                                        nPixelsArray);

            brainCRF{i}.imagesnode = [ones(1,1,size(images{i},3)) images{i} ...
                                UGM_standardizeCols(nBors{i},tied), ...
                                c1priors];
            sharedFeatures = [1 0 1 0];
        elseif (useCdist)
            % We know that useCdist and usePriors are
            % mututally exclusive
            brainCRF{i}.imagesnode = [ones(1,1,size(images{i},3)) images{i} ...
                                UGM_standardizeCols(nBors{i},tied) ...
                                cDist{i}];

            %add feature matricies
            sharedFeatures = [1 0 1 1];
            
        else
            brainCRF{i}.imagesnode = [ones(1,1,size(images{i},3)) images{i} ...
                                UGM_standardizeCols(nBors{i},tied)];

            %add feature matricies
            sharedFeatures = [1 0 1];
            
        end
        
        brainCRF{i}.imagesedge = ...
            UGM_makeEdgeFeatures(brainCRF{i}.imagesnode, ...
                                 brainCRF{i}.edgeStruct ...
                                 .edgeEnds,sharedFeatures(:));
        
        %Makes mapping of features to parameters
        [brainCRF{i}.nodeMap brainCRF{i}.edgeMap] = ...
            UGM_makeCRFmaps(brainCRF{i}.imagesnode, ...
                            brainCRF{i}.imagesedge, ...
                            brainCRF{i}.edgeStruct,0,tied,1,1);
        
        brainCRF{i} = save_data(dir, brainCRF{i}, i);
        
        if (usePriors)
            clear('c1priors')
            %clear('c2priors')
            %clear('c3priors')
        end
    end
    fprintf('\n');
end


function priors = load_spm8_matrix(res, tissueNum, imageNum, Zmask, ...
                                   nPixelsArray)
    
    priorsx = load_spm8_priors(res, tissueNum,imageNum);
    
    restarting = 0;
    
    nPixels = nPixelsArray(imageNum);
    
    if restarting == 0
        
        priorsx = priorsx(Zmask{imageNum});
        priorsx = reshape(priorsx,1,1,nPixels);
    end
    
    priors = priorsx;
end

function images = load_spm8_priors(res, tissueNum, imageNum)
% We will only use one of either images or y, but we are unsure of which
% to use as of yet    
% IF WE WANT TO USE THIS, WE HAVE TO CHANGE THIS DIRECTORIES
    file = '/acmi/fmri/IBSR_nifti_stripped/new_Segment_MRF2_dist2/';
    % Directory we use for the
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
    images = I_T1;
    

    %fprintf('%d',size(images,3),size(y,3))
    images = images(1:res:end,1:res:end,1:res:end);
    
end

function decode(weights, brainCRF, ZmaskFlat, origimages, numImages, indImages)
    
    
    for i = 1:numImages

        brainCRF{i} = load_data(brainCRF{i});
        [nodePot,edgePot] = UGM_CRF_makePotentials(weights, ...
                                                   brainCRF{i}.imagesnode, ...
                                                   brainCRF{i}.imagesedge, ...
                                                   brainCRF{i}.nodeMap, ...
                                                   brainCRF{i}.edgeMap, ...
                                                   brainCRF{i}.edgeStruct,1);
        
        yDecode = UGM_Decode_ICMrestart(nodePot,edgePot,...
                                        brainCRF{i}.edgeStruct,30); %last value is number of restarts
        disp('Size of ydecode:')
        disp(size(yDecode));
        disp('Size of ZmaskFlat{i}:');
        disp(size(ZmaskFlat{i}))
        yDecode = reImage(yDecode, ZmaskFlat{i});
        disp('Size of ydecode:')
        disp(size(yDecode));        
        yDecode(yDecode == 0) = 1;
        [nRows, nCols, nSlices] = size(origimages{i});
        disp('Size of origimages{i}:');
        disp(size(origimages{i}));
        yDecode = reshape(yDecode, nRows, nCols, nSlices);
        disp('Size of ydecode:')
        disp(size(yDecode));         

        imOut = make_nii(yDecode);
        save_nii(imOut, strcat('/acmi/summer2014/ADNI_tissues/', ...
                               type,sprintf('%03d',indImages(i)), ...
                               '_tissueSeg.nii'));
        brainCRF{i} = save_data(dir, brainCRF{i}, i);
    end
end


function cDist = center_distance(images,numImages)
    cDist = cell(numImages,1);
    for i=1:numImages
        xlen = size(images{i},1);
        ylen = size(images{i},2);
        zlen = size(images{i},3);
        [A,B,C] = meshgrid(1:xlen,1:ylen,1:zlen);
        cDist{i}=sqrt(((A-xlen/2).^2)+((B-ylen/2).^2)+((C-zlen).^2));
    end
    cDist = cor_bias(cDist,numImages);
end

function nImage = reImage(masked, ZmaskFlat)
    nImage = int32(ZmaskFlat);
    nImage(:,:,ZmaskFlat) = masked(:,:,:);
end
