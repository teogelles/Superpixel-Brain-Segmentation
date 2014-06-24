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
    centerAddr=strcat(saveDir,'centers','-',dirType,'-', ...
                      int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                      '-',int2str(res),'-',int2str(imageNum), ...
                      '-',int2str(numIters),'.mat');
    trackerAddr=strcat(saveDir,'tracker','-',dirType,'-', ...
                       int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                       '-',int2str(res),'-',int2str(imageNum), ...
                       '-',int2str(numIters),'.mat');
    indexAddr=strcat(saveDir,'index','-',dirType,'-', ...
                       int2str(numSuperVoxels),'-',int2str(shapeParam), ...
                       '-',int2str(res),'-',int2str(imageNum), ...
                       '-',int2str(numIters),'.mat');
    
    
    % checks if we've already run our primary SLIC code and thus
    % the file already exists
    if (exist(slicAddr, 'file') && exist(borderAddr, 'file') && ...
        exist(xAddr, 'file') && exist(centerAddr, 'file') && ...
        exist(trackerAddr, 'file') && exist(indexAddr, 'file'))
        
        fprintf('Relevant Files Already Exist, Loading...\n');
        
        labels = load_nifti(slicAddr, imageNum, 1);
        centers = load(centerAddr);
        centerTracker = load(trackerAddr);
        indexList = load(indexAddr);
        
        centers = centers.centers;
        centerTracker = centerTracker.centerTracker;
        indexList = indexList.indexList;
        
        tissueFilename = strcat('/acmi/chris13/results/ADNIresults/', ...
                                dirType, int2str(imageNum), '_again');

        featureFilename = '/';
        slicFeatures = getSLICFeatures(labels, centers, centerTracker, ...
                                               tissueFilename, ...
                                               indexList, featureFilename);
    else
        
        [X indexList] = load_nifti(dirType,imageNum,res);
        
        [labels border centers centerTracker] = SLIC_3D(X,numSuperVoxels, ...
                                                        shapeParam, numIters);
        
        slicNii = make_nii(labels);
        borderNii = make_nii(border);
        xNii = make_nii(X);
        
        fprintf('Saving SLIC To %s\n', slicAddr);
        fprintf('Saving Border to %s\n', borderAddr);
        fprintf('Saving X to %s\n', xAddr);
        fprintf('Saving Centers to %s\n', centerAddr);
        fprintf('Saving Tracker to %s\n', trackerAddr);
        
        save_nii(slicNii, slicAddr);
        save_nii(borderNii, borderAddr);
        save_nii(xNii, xAddr);
        save(centerAddr, 'centers');
        save(trackerAddr, 'centerTracker');
        save(indexAddr, 'indexList');
        
        tissueFilename = strcat('/acmi/chris13/results/ADNIresults/', ...
                                dirType, int2str(imageNum), ...
                                '_again');
        
        
        featureFilename = '/';
        slicFeatures = getSLICFeatures(labels, centers, centerTracker, ...
                                               tissueFilename, ...
                                               indexList, featureFilename, ...
                                               res);
    end
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
    
    
    
    indexList = [1, size(X, 1); 1, size(X, 2); 1, size(X, 3)];
    [X indexList] = cropBlack(X);
end