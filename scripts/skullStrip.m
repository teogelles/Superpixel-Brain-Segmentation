% Attempt at skull stripping program for flain and MRI image data

function skullStrip(threshold, filename, outfilename)
   if ~exist('filename', 'var')
       filename = '/sonigroup/chris13/data/DHS/data/16902.nii';
   end
   
   if ~exist('outfilename','var')
       outfilename = 'stripped.nii';
   end
   
   if ~exist('threshold','var')
       % Must be between 0 and 1
       threshold = .2;
   end
   
   mri = load_nifti(filename);
   minIntensity = min(mri(:));
   maxIntensity = max(mri(:));
   % tIn = thresholdIntensity
   tIn = minIntensity + (maxIntensity - minIntensity)*threshold;
    
    
    
    
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
