function normalizeImage(filename, outfilename)
    
    X = load_nifti(filename);
 
    %X = load_nii(filename);
    
    length = size(X, 1);
    width = size(X, 2);
    pages = size(X, 3);
    
    for i=1:length
        for j=1:width;
            for k=1:pages
                
                X(i, j, k) = X(i, j, k) + .1*X(i, j, k);
            end
        end
    end
    
    outfilename = '/acmi/summer2014/tgelles1/brainseg2014/scripts/blank.nii';
    save_nii(X, outfilename);
end

function X = load_nifti(filename)
%Loads IBSR V2 nifti files

    fprintf('Loading Nifti Image\n');
    
    %Image
    I_t1uncompress = wfu_uncompress_nifti(filename);
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X = int32(I_T1);
end