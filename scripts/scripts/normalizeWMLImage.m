function normalizeWMLImage(imNum, outfilename)
    
    
    if ~exist('imNum','var')
        imNum = 16902;
    end
    
    imNum = num2str(imNum);
    
    if ~exist('outfilename','var')
        outfilename = strcat('/sonigroup/summer2014/agilchr1/brainseg2014/', ...
                       'scripts/stripped_', imNum, '.nii');
    end
        
    filename = strcat('/sonigroup/chris13/data/DHS/test/',imNum, ...
                      '.nii');
    c1filename = strcat('/sonigroup/chris13/data/DHS/test/c1',imNum, ...
                      '.nii');
    c2filename = strcat('/sonigroup/chris13/data/DHS/test/c2',imNum, ...
                      '.nii');
    c3filename = strcat('/sonigroup/chris13/data/DHS/test/c3',imNum, ...
                      '.nii');
        
    X = load_nifti(filename);
    c1 = load_nifti(c1filename);
    c2 = load_nifti(c2filename);
    c3 = load_nifti(c3filename);
    
    thresh = .1;
    X = skullStrip(X,c1,c2,c3,thresh);
    
    %X = load_nii(filename);
    
    % length = size(X, 1);
    % width = size(X, 2);
    % pages = size(X, 3);
    
    % for i=1:length
    %     for j=1:width;
    %         for k=1:pages
    %             nebEnds = getNeighborhoodEnds(size(X),1,i,j,k);
    %             wiper = true;
    %             for ii = nebEnds(1):nebEnds(2)
    %                 for jj = nebEnds(3):nebEnds(4)
    %                     for kk = nebEnds(5):nebEnds(6)
    %                         if X(i,j,k) ~= 0
    %                             wiper = false;
    %                         end
    %                     end
    %                 end
    %             end
    %             if wiper
    %                 X(i,j,k) = 0;
    %             end
    %         end
    %     end
    % end
    
    Xnii = make_nii(X);
    fprintf('Saving output to %s\n', outfilename);
    save_nii(Xnii, outfilename);
    
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


function neighborhoodEnds = getNeighborhoodEnds(imageMatSize, radius, i, j, k)
% Function gets the neighborhood ends for a regions around the
% center of size radius
% @param imageMatSize - size of the image matrix
% @param radius - square radius of neighborhood around center
% @params i,j,k - x,y,z coordintate of center
%
% @return neighborhoodEnds - matrix of starts and ends of each
% direction of neighborhood
    
    neighborhoodEnds = [floor(i-radius),ceil(i+radius),floor(j-radius), ...
                        ceil(j + radius), floor(k - radius), ceil(k + radius)];
    
    % edge cases
    if neighborhoodEnds(1) < 1
        neighborhoodEnds(1) = 1;
    end
    if neighborhoodEnds(3) < 1
        neighborhoodEnds(3) = 1;
    end
    if neighborhoodEnds(5) < 1
        neighborhoodEnds(5) = 1;
    end
    
    if neighborhoodEnds(2) > imageMatSize(1)
        neighborhoodEnds(2) = imageMatSize(1);
    end
    if neighborhoodEnds(4) > imageMatSize(2);
        neighborhoodEnds(4) = imageMatSize(2);
    end
    if neighborhoodEnds(6) > imageMatSize(3);
        neighborhoodEnds(6) = imageMatSize(3);
    end
    
end

function X = skullStrip(X,c1,c2,c3,threshold);
    fprintf('Masking skull\n');
    mask = ((c1 + c2 + c3) > threshold);
    mask = int32(mask);
    X = X.*mask;
end
