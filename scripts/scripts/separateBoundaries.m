%Attempt at a script which takes our boundary images and removes
%the boundaries from the original image background

function separateBoundaries(type,imNum)
    
    if ~exist('type','var')
        type = 'AD';
    end
    
    if ~exist('imNum','var')
        imNum = 1;
    end
    
    filebase = '/scratch/tgelles1/summer2014/slic/';
    
    borderAddr = strcat('border-',type,'-500-0.1-1-', ...
                        num2str(imNum),'-18.nii');
    fullborderAddr = strcat(filebase,borderAddr);
    origAddr =  strcat('x-500-0.1-',type,'-1-', ...
                        num2str(imNum),'-18.nii');
    fullorigAddr = strcat(filebase,origAddr);
    
    borders = load_nifti(fullborderAddr);
    orig = load_nifti(fullorigAddr);
    
    justBorderInd = find(~borders & orig);
    justBorder = zeros(size(borders));
    justBorder(justBorderInd) = 1;
    
    justB_nii = make_nii(justBorder);
    justBname = strcat(filebase,'just',borderAddr);
    save_nii(justB_nii,justBname);
    
end

function X = load_nifti(imageName);
    
    I_t1uncompress = wfu_uncompress_nifti(imageName);
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X = I_T1;
end