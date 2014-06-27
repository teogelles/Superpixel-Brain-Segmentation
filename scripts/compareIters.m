% simple file to compare the effect of the number of iterations has
% on runSLIC

function compareIters()
    
    base = '/scratch/tgelles1/summer2014/slic_test/slic-IBSR-500-20-1-1-';
    
    %Set for first iteration, then old will get overwritten, making
    %the name make more sense
    oldLabels = load_nifti(strcat(base,'4.nii'));
    
    numVox = size(oldLabels,1)*size(oldLabels,2)*size(oldLabels,3);
    
    for i = 5:30
        filename = strcat(base,num2str(i),'.nii');
        newLabels = load_nifti(filename);
        E = sum(sum(sum(newLabels ~= oldLabels)));
        percentE = (E/numVox);
        fprintf('%d change between %d and %d iterations\n', ...
                percentE,i,i-1);
        oldLabels = newLabels;
    end
end

function ret = load_nifti(filename)
    I_t1uncompress = wfu_uncompress_nifti(filename);
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    ret = I_T1;
end