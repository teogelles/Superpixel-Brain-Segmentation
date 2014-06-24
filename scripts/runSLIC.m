function slicFeatures = runSLIC(imageNum, dirType, res, numSuperVoxels, ...
                                shapeParam, numIters)
    
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
        
        fprintf('Saving SLIC To %s\n', slicAddr);
        fprintf('Saving Border to %s\n', borderAddr);
        fprintf('Saving X to %s\n', xAddr);
        fprintf('Saving CenterInfo to %s\n', centerinfoAddr);
        
        save_nii(slicNii, slicAddr);
        save_nii(borderNii, borderAddr);
        save_nii(xNii, xAddr);
        save(centerinfoAddr, 'centerInfo');
        save(cropAddr, 'cropOffset');
    end
    
    
    if (strcmp(dirType, 'IBSR'))
        tissues = NaN;

    else
        tissueFilename = strcat('/acmi/chris13/results/ADNIresults/', ...
                                dirType, int2str(imageNum), '_again');
        tissues = load_tissues(tissueFilename, cropOffset, res);
    end
    
    
    featureFilename = '/';
    slicFeatures = getSLICFeatures(labels, tissues, centerInfo, ...
                                           cropOffset, featureFilename);
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
    elseif (strcmp(dirType, 'MCN'))
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