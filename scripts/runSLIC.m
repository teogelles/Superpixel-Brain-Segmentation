% slicFeatures
% Authors: Andrew Gilchrist-Scott & Teo Gelles
%
% This file contains the code for runSLIC, which is the main
% wrapper used in MATLAB in order to run experiments with SLIC_3D
% and getSLICFeatures.

function slicFeatures = runSLIC(imageNum, dirType, res, numSuperVoxels, ...
                                shapeParam, numIters)
    % slicFeatures - Returns the list of features obtained from
    % getSLICFeatures()
    %
    % @param imageNum - The numerical index of the image to use in
    % its given folder
    % @param dirType - The basic directory under which the image
    % can be found.  Currently recognizes 'CN', 'MCI', 'AD', or 'IBSR',
    % which refer to /acmi/fmri/CN_T1, /acmi/fmri/MCI_T1,
    % /acmi/fmri/AD_T1, and /acmi/fmri/IBSR_nifti_stripped
    % respectively.   Inputs which do not match one of these are
    % assumed to refer to the entire directory for the given image.
    % @param res - The inverse resolution of the image (1 for full,
    % 2 for half, etc.)
    % @param numSuperVoxels - The number of superVoxels to use in
    % SLIC
    % @param shapeParam - The weight to use for the distance metric
    % in SLIC
    % @param numIters - The number of iterations to run SLIC for
    
    % Handles if the user chooses to not input any of the arguments
    if ~exist('numIters','var')
        numIters = 15;
    end
    
    if ~exist('shapParam','var')
        shapeParam = 20;
    end
    
    if ~exist('numSuperVoxels','var')
        numSuperVoxels = 500;
    end
    
    if ~exist('res','var')
        res = 1;
    end
    
    if ~exist('dirType','var')
        dirType = 'IBSR';
    end
    
    if ~exist('imageNum','var')
        imageNum = 1;
    end
    
    % base directory
    saveDir = '/scratch/tgelles1/summer2014/slic/';
    
    % file addressing specific to each of the different type of
    % file we may choose to run
    slicAddr=strcat(saveDir,'slic','-',dirType,'-', ...
                    int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                    '-',int2str(res),'-',int2str(imageNum),'-', ...
                    int2str(numIters),'.nii');
    borderAddr=strcat(saveDir,'border','-',dirType,'-', ...
                      int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                      '-',int2str(res),'-',int2str(imageNum), ...
                      '-',int2str(numIters),'.nii');
    xAddr=strcat(saveDir,'x','-',int2str(numSuperVoxels),'-', ...
                 int2str(shapeParam),'-',dirType,'-',int2str(res), ...
                 '-',int2str(imageNum),'-',int2str(numIters),'.nii');
    centerinfoAddr=strcat(saveDir,'centerinfo','-',dirType,'-', ...
                          int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                          '-',int2str(res),'-',int2str(imageNum), ...
                          '-',int2str(numIters),'.mat');
    cropAddr=strcat(saveDir,'cropoffset','-',dirType,'-', ...
                    int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                    '-',int2str(res),'-',int2str(imageNum), ...
                    '-',int2str(numIters),'.mat');
    
    
    % checks if we've already run our primary SLIC code and thus
    % the file already exists
    if (exist(slicAddr, 'file') && exist(borderAddr, 'file') && ...
        exist(xAddr, 'file') && exist(centerinfoAddr, 'file') && ...
        exist(cropAddr, 'file'))
        
        fprintf('Relevant Files Already Exist, Loading...\n');
        
        labels = load_nifti(slicAddr, imageNum, 1);
        X = load_nifti(xAddr, imageNum, 1);
        
        centerInfo = load(centerinfoAddr);
        cropOffset = load(cropAddr);
        
        centerInfo = centerInfo.centerInfo;
        cropOffset = cropOffset.cropOffset;
    else
        
        [X cropOffset] = load_nifti(dirType,imageNum,res);
        
        [labels border centerInfo] = SLIC_3D(X,numSuperVoxels, ...
                                             shapeParam, numIters);
        
        slicNii = make_nii(labels);
        borderNii = make_nii(border);
        xNii = make_nii(X);
        
        fprintf('Saving SLIC to %s\n', slicAddr);
        save_nii(slicNii, slicAddr);
        fprintf('Saving Border to %s\n', borderAddr);
        save_nii(borderNii, borderAddr);
        fprintf('Saving X to %s\n', xAddr);
        save_nii(xNii, xAddr);
        fprintf('Saving CenterInfo to %s\n', centerinfoAddr);
        save(centerinfoAddr, 'centerInfo');
        fprintf('Saving CropOffset to %s\n', cropAddr);
        save(cropAddr, 'cropOffset');
    end
    
    
    if (strcmp(dirType, 'IBSR'))
        tissues = NaN;

    else
        tissueFilename = strcat('/acmi/chris13/results/ADNIresults/', ...
                                dirType, int2str(imageNum), '_again');
        tissues = load_tissues(tissueFilename, cropOffset, res);
    end
    
    if strcmp(dirType, 'AD')
        id = imageNum;
    elseif strcmp(dirType, 'MCI')
        id = imageNum + 92;
    elseif strcmp(dirType, 'CN')
        id = imageNum + 92 + 203;
    else
        id = imageNum+1000;
    end
    
    featureFilename = strcat('/scratch/tgelles1/summer2014/slic/',...
                             dirType,num2str(imageNum),'.txt');
    
    slicFeatures = getSLICFeatures(X, labels, tissues, centerInfo, ...
                                      cropOffset,featureFilename, ...
                                      id);
end


function [X, indexList] = load_nifti(dirType,imageNum, res)


    fprintf('Loading Nifti Image...\n');
    
    imageName = '';
    if (strcmp(dirType, 'IBSR'))
        
        if (imageNum < 10)
            imageName = strcat('/acmi/fmri/IBSR_nifti_stripped/IBSR_0', ...
                               int2str(imageNum), '_ana_strip.nii');
        else
            imageName = strcat('/acmi/fmri/IBSR_nifti_stripped/IBSR_', ...
                               int2str(imageNum), '_ana_strip.nii');
        end
        
    elseif (strcmp(dirType, 'AD'))
        
        imageName = strcat('/acmi/fmri/AD_T1/patient', ...
                           int2str(imageNum), '.nii');
    elseif (strcmp(dirType, 'CN'))
        imageName = strcat('/acmi/fmri/CN_T1/patient', ...
                           int2str(imageNum), '.nii');
    elseif (strcmp(dirType, 'MCI'))
        imageName = strcat('/acmi/fmri/MCI_T1/patient', ...
                           int2str(imageNum), '.nii');
    else
        imageName = dirType;
    end

    if (~exist(imageName, 'file'))
        exception = MException('file:dne', ['file %s does not exist ' ...
                            'or you do not have proper permissions'], ...
                               imageName);
        throw(exception);
    end
    

    fprintf('ImageName: %s\n', imageName);
    
    %Image
    I_t1uncompress = wfu_uncompress_nifti(imageName);
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X = I_T1;
    
    
    X = X(1:res:end,1:res:end,1:res:end);
    [X indexList] = cropBlack(X);
end

function tissues = load_tissues(tissueFilename, cropOffset, res)
    
    tissues = load_nii(tissueFilename);
    tissues = tissues.img(1:res:end,1:res:end,1:res:end);
    tissues = tissues(cropOffset(1, 1):cropOffset(1, 2), cropOffset(2, ...
                                                      1):cropOffset(2, ...
                                                      2), cropOffset(3, ...
                                                      1):cropOffset(3, ...
                                                      2));
end