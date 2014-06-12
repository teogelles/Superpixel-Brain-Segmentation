function convertNifti()
[X, y, nExamples] = load_ibsr('/acmi/chris13/data/analyze_IB/20Normals_T1_brain/', '/acmi/chris13/data/20Normals_T1_seg/'); %Load IBSR_V1
fprintf('Saving files');
for i=1:nExamples
    fprintf('.');
    X{i} = permute(X{i}, [1 3 2]);
    y{i} = permute(y{i}, [1 3 2]);
    %X{i} = flipdim(X{i},1);
    %X{i} = flipdim(X{i},2);
    X{i} = flipdim(X{i},3);
    %y{i} = flipdim(y{i},1);
    %y{i} = flipdim(y{i},2);
    y{i} = flipdim(y{i},3);
    tempNiiX = make_nii(X{i});
    tempNiiY = make_nii(y{i});
    save_nii(tempNiiX, strcat('/acmi/chris13/data/ibsr_v1_nii/im',int2str(i)));
    save_nii(tempNiiY, strcat('/acmi/chris13/data/ibsr_v1_nii/seg',int2str(i)));
end
fprintf('Success!\n');
end


function [X,y,nExamples] = load_ibsr(imDir, segDir)
offsets = [0 1 2 1 1 0 0 0 0 3 3 -4 3 2 6 8 1 0 0 2];
%Load data from ibsr database
bList = dir(strcat(imDir,'*.hdr'));
rawImages = cell(length(bList));
nExamples = length(bList)-1; %CAUSE WHAT THE HELL -4



bList = dir(strcat(segDir,'*.hdr'));
segs = cell(length(bList));
ends = zeros(nExamples,1);

for f = 1:length(bList)
    fid = fopen(strcat(segDir,bList(f).name));
    stats = fscanf(fid, '%d');
    nSlices = stats(3);
    ends(f) = nSlices;
    fid = fopen(strcat(segDir,bList(f).name(1:end-4),'.buchar'));
    I = fread(fid,inf,'uint8=>int32');
    I = reshape(I, 256, 256, nSlices);
    if offsets(f) > -1
        segs{f} = I;
    end
end

for f = 1:length(bList)
    fid = fopen(strcat(imDir,bList(f).name));
    stats = fscanf(fid, '%d');
    nSlices = stats(3)+2;
    fid = fopen(strcat(imDir,bList(f).name(1:end-4),'.buchar'));
    I = fread(fid,inf,'uint8=>int32');
    I = reshape(I, 256, 256, nSlices);
    if (ends(f) + offsets(f)) > nSlices
        fprintf('herr!, %d, %d\n', nSlices, ends(f));
        a = segs{f};
        segs{f} = a(:,:,1:nSlices-offsets(f));
        rawImages{f} = I(:,:,offsets(f)+1:end);
    else
        if offsets(f) > -1
            rawImages{f} = I(:,:,offsets(f)+1:ends(f)+offsets(f));
        end
    end
end

segs{12} = segs{11};
rawImages{12} = rawImages{11};

fclose('all'); %python needs this

%Create X and Y
y = segs;
X = rawImages;
for i=1:nExamples
    %1 is background, 2 is CSF, 3 is GM, 4 is WM
    y{i} = ((y{i}==0) + ((y{i}==128)*2) + ((y{i}==192)*3) + ((y{i}==254)*4)); 
    X{i} = double(X{i});
    y{i} = int32(y{i});
end
for i = 1:nExamples
    fprintf('%d, %d, %d\n',size(X{i},3),size(y{i},3), offsets(i))
    X{i} = X{i}(1:2:end,1:4:end,1:2:end);
     y{i} = y{i}(1:2:end,1:4:end,1:2:end);
end
% nExamples = 7;
% X(8:end) = [];
% y(8:end) = [];
end