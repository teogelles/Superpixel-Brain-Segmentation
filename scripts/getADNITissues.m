%% getADNITissues
% This program uses a CRF to determine the tissue segmentation on
% ~300 ADNI brain images.
%
% Written by Teo Gelles and Andrew Gilchrist-Scott
% Advised by Ameet Soni

function getADNITissues(type,section,wegihtFile,useCdist,usePriors)
% @param type - Type of subject (AD,MCI, or CN)
% @param section - portion of each types data we want to run. This
% is useful for using condor on this functoin sinc we can decode
% multiple images simultaneously. Should be in [1,10]
% @param weightFile - filename of .mat file with the weights from
% the optimal CRF training using the full 18 IBSR images
    
    checkInput(type,section,weightfile);
    
    global dir;
    dir = makeDirs(usePriors, useCdist);
        
    % Load the data
    weights = load(weightfile);
    weights = weights.w;
    images = load_nifti(type,section);
    numImages = size(images,1);
    
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
    [origimages, origY, Zmask, ZmaskFlat, images, y, nStates, sizes, nPixelsArray, ...
     nBors, cDist] = reshapeMatrices(numImages, images, nBors, cDist, ...
                                     useCdist);

    
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

function images = load_nifti(type,section)
%Loads IBSR V2 nifti files

    fprintf('Loading Nifti Images');
    
    patients = getPatients(type,section);
    
    images = cell(size(patients,1),1);
    
    parfor pat_i = pat_iients
        
        fprintf('.');
        filename = strcat(['/scratch/tgelles1/summer2014/' ...
                           'ADNI_stripped'],type,pat_i,'.nii');
        
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
            patientsOfType = totalPatients{i};
            break
        end
    end
    
    step = floor(patientsOfType/9);
    if section == 1
        patients = 1:step
    elseif section == 10
        patients = (step*9):patientsOfType;
    else
        patients = (section*(step-1) + 1):(section*step);
    end
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

function [origimages, origY, Zmask, ZmaskFlat, images, y, nStates, sizes, ...
          nPixelsArray, nBors, cDist] = reshapeMatrices(numImages, images, nBors, cDist, ...
                                     useCdist);
    
    global dir;
    
    sizes = zeros(numImages);
    parfor i=1:numImages
        [nRows,nCols,nSlices] = size(images{i});
        sizes(i) = nRows*nCols*nSlices;
    end
    
    %%%STOPPED WORKING HERE 6/27a

    fprintf('Masking Zeros...\n');
    [origimages, origY, Zmask, images, y] = maskZeros(images, y, numImages);

    nStates = max(y{1}(:)); %assume y{1} has all states 
    ZmaskFlat = cell(numImages,1);
    
    nPixelsArray = zeros(1, numImages);

    fprintf('Reshaping Matricies');
    for i=1:numImages
        fprintf('.');
        nPixels = size(images{i},1);
        nPixelsArray(i) = nPixels;
        
        y{i} = reshape(y{i},[1, 1 nPixels]);
        images{i} = reshape(images{i},1,1,nPixels);
        
        if restarting == 0
            nBors{i} = nBors{i}(Zmask{i});
            nBors{i} = reshape(nBors{i},1,1,nPixels);
            
            if (useCdist)
                cDist{i} = cDist{i}(Zmask{i});
                cDist{i} = reshape(cDist{i},1,1,nPixels);
            end
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