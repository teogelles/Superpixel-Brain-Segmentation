function slicFeatures = runSLIC(res, numSuperVoxels, shapeParam, ...
                                numIters, dirType, imageNum)
    
    saveDir = '/scratch/tgelles1/summer2014/slic/';
    
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
    centerAddr=strcat(saveDir,'centers','-',dirType,'-', ...
                      int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                      '-',int2str(res),'-',int2str(imageNum), ...
                      '-',int2str(numIters),'.mat');
    trackerAddr=strcat(saveDir,'tracker','-',dirType,'-', ...
                       int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                       '-',int2str(res),'-',int2str(imageNum), ...
                       '-',int2str(numIters),'.mat');
    
    
    
    if (exist(slicAddr, 'file') && exist(borderAddr, 'file') && ...
        exist(xAddr, 'file') && exist(centerAddr, 'file') && ...
        exist(trackerAddr, 'file'))
        
        
        fprintf('Relevant Files Already Exist, Loading...\n');
        
        labels = load_nifti(slicAddr, imageNum, 1);
        centers = load(centerAddr);
        centerTracker = load(trackerAddr);
        
        centers = centers.centers;
        centerTracker = centerTracker.centerTracker;
        
        featureFilename = '/';
        slicFeatures = getSLICFeatures(labels, centers, centerTracker, ...
                                               featureFilename);
    else
        
        X = load_nifti(dirType,imageNum,res);
        
        [labels border centers centerTracker] = SLIC_3D(X,numSuperVoxels, ...
                                                        shapeParam, numIters);
        
        slicNii = make_nii(labels);
        borderNii = make_nii(border);
        xNii = make_nii(X);
        
        fprintf('Saving SLIC To %s\n', slicAddr);
        fprintf('Saving Border to %s\n', borderAddr);
        fprintf('Saving X to %s\n', xAddr);
        fprintf('Saving Centesr to %s\n', centerAddr);
        fprintf('Saving Tracker to %s\n', trackerAddr);
        
        save_nii(slicNii, slicAddr);
        save_nii(borderNii, borderAddr);
        save_nii(xNii, xAddr);
        save(centerAddr, 'centers');
        save(trackerAddr, 'centerTracker');
        
        featureFilename = '/';
        slicFeatures = getSLICFeatures(labels, centers, centerTracker, ...
                                               featureFilename);
    end
end


function X = load_nifti(dirType,imageNum, res)


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
    X = cropBlack(X);
end