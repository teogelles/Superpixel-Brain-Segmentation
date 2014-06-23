function runSLIC(res, numSuperVoxels, shapeParam, numIters, imageNum)
    
    
    X = load_nifti('/acmi/fmri/IBSR_nifti_stripped/', ...
                   res, imageNum);
    
    [labels border centers centerTracker] = SLIC_3D(X,numSuperVoxels, ...
                                                      shapeParam, numIters);
        
    slicNii = make_nii(labels);
    borderNii = make_nii(border);
    xNii = make_nii(X);
    
    saveDir = '/scratch/tgelles1/summer2014/slic/';
    
    slicAddr = strcat(saveDir, 'slic', '-', int2str(numSuperVoxels), ...
                             '-', int2str(shapeParam), '-', int2str(res), ...
                             '-', int2str(imageNum), '-', int2str(numIters), '.nii');
    borderAddr = strcat(saveDir, 'border', '-', int2str(numSuperVoxels), ...
                             '-', int2str(shapeParam), '-', int2str(res), ...
                             '-', int2str(imageNum), '-', int2str(numIters), '.nii');
    xAddr = strcat(saveDir, 'x', '-', int2str(numSuperVoxels), ...
                             '-', int2str(shapeParam), '-', int2str(res), ...
                             '-', int2str(imageNum), '-', int2str(numIters),'.nii');
    
    fprintf('Saving SLIC To %s\n', slicAddr);
    fprintf('Saving Border to %s\n', borderAddr);
    fprintf('Saving X to %s\n', xAddr);
    
    save_nii(slicNii, slicAddr);
    save_nii(borderNii, borderAddr);
    save_nii(xNii, xAddr);
    
    
    featureFilename = '/';
    getSLICFeatures(labels, centers, centerTracker, featureFilename);
end


function X = load_nifti(imDir,res, imageNum)


    fprintf('Loading Nifti Image...\n');
    
    if imageNum < 10
        place = strcat(imDir,'IBSR_0',int2str(imageNum),'/');
        fileHead = strcat('IBSR_0',int2str(imageNum));
    else
        place = strcat(imDir,'IBSR_',int2str(imageNum),'/');
        fileHead = strcat('IBSR_',int2str(imageNum));
    end
    
    %Image
    I_t1uncompress = wfu_uncompress_nifti(strcat(place,fileHead,'_ana_strip.nii'));
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X = I_T1;
    
    
    X = X(1:res:end,1:res:end,1:res:end);
    X = cropBlack(X);
end