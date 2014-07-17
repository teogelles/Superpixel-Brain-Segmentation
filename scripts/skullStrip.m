% Attempt at skull stripping program for flair and MRI image data

function skullStrip(threshold, filename, outfilename)
   if ~exist('filename', 'var')
       filename = '/sonigroup/chris13/data/DHS/data/16902.nii';
   end
   
   if ~exist('outfilename','var')
       outfilename = strcat('stripped',num2str(threshold),'.nii');
   end
   
   if ~exist('threshold','var')
       % Must be between 0 and 1
       threshold = .2;
   end
   
   mri = load_nifti(filename);
   mri = mri*2;
   minIntensity = min(mri(:));
   maxIntensity = max(mri(:));
   % tIn = thresholdIntensity
   tIn = minIntensity + (maxIntensity - minIntensity)*threshold;
   mri = strip(mri, tIn);
   mriNii = make_nii(mri);
   save_nii(mriNii, outfilename);
   
   fprintf('Done\n');
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

function mri = strip(mri, tIn) 
    
    fprintf('Stripping');
    
    inSkull = false;
    outSkull = false;
    headMask = uint8(mri > tIn);
    notHeadMask = ~headMask;
    for i = 1:size(mri,1)
        if ~mod(i,10)
            fprintf('.');
        end
        for j = 1:size(mri,2)
            lit = find(headMask(i,j,:));
            if size(lit(:),1) == 0
                % line in matrix doesn't include head
                mri(i,j,:) = 0;
                continue
            end
            darkInSkull = find(notHeadMask(i,j,lit(1):lit(end)));
            if size(darkInSkull(:),1) == 0
                % line in matrix only includes skull
                mri(i,j,:) = 0;
                continue
            end
            mri(i,j,1:darkInSkull(1)) = 0;
            mri(i,j,darkInSkull(end):end) = 0;
        end
    end
    fprintf('\n'); 
    for i = 1:size(mri,1)
        if ~mod(i,10)
            fprintf('.');
        end
        for k = 1:size(mri,3)
            lit = find(headMask(i,:,k));
            if size(lit(:),1) == 0
                % line in matrix doesn't include head
                mri(i,:,k) = 0;
                continue
            end
            darkInSkull = find(notHeadMask(i,lit(1):lit(end),k));
            if size(darkInSkull(:),1) == 0
                % line in matrix only includes skull
                mri(i,:,k) = 0;
                continue
            end
            mri(i,1:darkInSkull(1),k) = 0;
            mri(i,darkInSkull(end):end,k) = 0;
        end
    end
    fprintf('\n');  
        for j = 1:size(mri,1)
        if ~mod(j,10)
            fprintf('.');
        end
        for k = 1:size(mri,3)
            lit = find(headMask(:,j,k));
            if size(lit(:),1) == 0
                % line in matrix doesn't include head
                mri(:,j,k) = 0;
                continue
            end
            darkInSkull = find(notHeadMask(lit(1):lit(end),j,k));
            if size(darkInSkull(:),1) == 0
                % line in matrix only includes skull
                mri(:,j,k) = 0;
                continue
            end
            mri(1:darkInSkull(1),j,k) = 0;
            mri(darkInSkull(end):end,j,k) = 0;
        end
    end
    fprintf('\n');  
end
