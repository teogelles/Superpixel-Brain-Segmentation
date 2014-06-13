function eval_seg()
aGM = 0;
aWM = 0;
aCSF = 0;
tWM = 0;
tGM= 0;
tCSF = 0;


%[segs,nSegs] = load_3File('/acmi/fmri/IBSR_nifti_stripped/coReg2/VBM/');
%[segs,nSegs] = load_3File('/acmi/chris13/data/ibsr_v1_nii/');
%[manual, nExamples] = load_1File('/acmi/chris13/v1_analyze/20Normals_T1_analyze/20Normals_T1_8bit/');
%[segs, nSegs] = load_fast_file('/acmi/chris13/v1_analyze/20Normals_T1_analyze/20Normals_T1_8bit/');
%[segs,nSegs] = load_3File('/acmi/fmri/IBSR_nifti_stripped/coReg2/oppss/');
%[segs, nSegs] = load_fast_file('/acmi/chris13/data/analyze_IB/FAST/');
%[nope, manual, nExamples] = load_ibsr_segs('/acmi/chris13/data/analyze_IB/20Normals_T1_brain/', '/acmi/chris13/v1_analyze/20Normals_T1_analyze/20Normals_T1_seg/'); %Load IBSR_V1
[manual, nExamples] = load_fast_file('/acmi/chris13/v1_analyze/20Normals_T1_analyze/20Normals_T1_seg/');
[segs,nSegs] = load_3File('/acmi/chris13/v1_analyze/20Normals_T1_analyze/');
%[manual,nExamples] = load_1File('/acmi/fmri/IBSR_nifti_stripped/');
%[bWeb, bWebMan, la] = load_brainweb();
%nExamples = 18;
%nSegs = 18;
segs{8} = segs{7};
manual{8} = manual{7};
% segs{1} = bWeb;
% manual{1} = bWebMan;
plots = 1;
if plots == 1
  X = manual;
      figure;
      lenT = nExamples;
      lenT = sqrt(lenT);
      for j = 1:nExamples
          subplot(lenT, lenT+1, j);
          imagesc(reshape((X{j}(:,:,floor(size(X{j},3)/2))),size(X{j},1),size(X{j},2))); %Should not hardcode these
          colormap gray
      end
      suptitle('MAN Images');
      
      
      X = segs;
      figure;
      lenT = nExamples;
      lenT = sqrt(lenT);
      for j = 1:nExamples
          subplot(lenT, lenT+1, j);
          imagesc(reshape((X{j}(:,:,floor(size(X{j},3)/2))),size(X{j},1),size(X{j},2))); %Should not hardcode these
          colormap gray
      end
      suptitle('MRI Images');
end
%offsets = [0 1 2 1 1 0 0 0 0 3 3 3 3 2 6 8 1 0 0 2];
%segs{12} =segs{11};
for i=1:nExamples
%     nSlices = size(manual{i},3);
%     if nSlices+offsets(i) > size(segs{i},3)
%        offsets(i) = abs(size(manual{i},3)-size(segs{i},3));
%        segs{i} = segs{i}(1:end,1:end,offsets(i)+1:end);  
%     else
%        segs{i} = segs{i}(1:end,1:end,offsets(i)+1:(offsets(i)+nSlices));
%     end
    i
    %tempNii = make_nii(segs{i});
    %save_nii(tempNii, strcat('/acmi/chris13/results/vbmSeg',int2str(i)));
    size(segs{i})
    size(manual{i})
    [results, tani] =  vOverlap(segs{i}, manual{i});
    aWM = aWM + results(1);
    aGM = aGM + results(2);
    aCSF = aCSF + results(3);
    tWM = tWM + tani(1);
    tGM = tGM + tani(2);
    tCSF = tCSF + tani(3);   
end
 fprintf('\n Average Volume Overlap: %f %f %f\n Average Tanimoto %f %f %f\n', aWM/nExamples, ...
    aGM/nExamples, aCSF/nExamples, tWM/nExamples, tGM/nExamples, tCSF/nExamples);   
end

function M = mcr(L,R)
%L is the labels, R is the results, returns the misclassification rate
dif = int32(L)-int32(R);
K = dif(:,:,:) ~= 0;
M = sum(sum(sum(K)));
M = (M/numel(R(R~=1)))*100;
fprintf('\nMCR: %f', M);
end

function [results tani] = vOverlap(L,R)
%L is the labels, R is the results, returns the % volume overlap
results = [0 0 0];
tani = [0 0 0];
fprintf('\n\n');
max(max(max(R)))
for i=2:max(max(max(R)))
    V1 = sum(sum(sum(L == i)));
    V2 = sum(sum(sum(R == i)));
    dif = ismember(L,i)+ismember(R,i);
    K = dif(:,:,:) == 2;
    N = sum(sum(sum(K)));
    M =((N*2)/(abs(V1+V2)))*100;
    results(i-1) = M;
    tani(i-1) =(N/(V1+V2-N));
    fprintf('\n%% Volume Overlap: %f, Tanimoto: %f', M, N/(V1+V2-N));
end
end

function M = vComp(L,R)
N = 0;
V1 = 0;
V2 = 0;
for i=2:max(max(max(R)))
    V1 = sum(sum(sum(L == i)));
    V2 = sum(sum(sum(R == i)));
    K = abs(V1-V2);
    N = sum(sum(sum(K)));
    M = (N/(abs(V1+V2)/2))*100;
    fprintf('\n%% Volume Difference: %f', M);
end
M = (N/(abs(V1+V2)/2))*100;

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
%     X{i} = X{i}(1:2:end,1:4:end,1:2:end);
%      y{i} = y{i}(1:2:end,1:4:end,1:2:end);
end
% nExamples = 7;
% X(8:end) = [];
% y(8:end) = [];
end

function [X,nExamples] = load_3File(imDir)

%load WM
%load GM
%load CSF

nExamples = 20;
X = cell(nExamples,1);

fprintf('\nLoading %d files',nExamples);
for i = 1:nExamples
    if i == 8
      continue;
    end
    fprintf('.');
    heads = {'c1','c2','c3'};
    results = cell(4);
    
    bList = dir(strcat(imDir,'r*.hdr'));
    
    for j=1:3
        %heads{j} = strcat(imDir,'IBSR_', num2str(i,'%02d'),'/', heads{j}, 'IBSR_', num2str(i,'%02d'), '_ana.nii');
        %heads{j} = strcat(imDir, heads{j}, 'mcoreg2IBSR_', num2str(i,'%02d'), '_ana.nii');
        %heads{j} = strcat(imDir,'IBSR_', num2str(i,'%02d'), '_ana_strip_pve_',num2str(i-1),'.nii');
        heads{j} = strcat(imDir, heads{j}, bList(i).name);
        I_t1uncompress = wfu_uncompress_nifti(heads{j});
        I_uncompt1 = spm_vol(I_t1uncompress);
        I_T1 = spm_read_vols(I_uncompt1);
        results{j} = int32(I_T1); 
    end
origSize = size(results{1});
marginals = zeros([3 size(results{2})]);
marginals(2,:,:,:) = results{1}(:,:,:);
marginals(3,:,:,:) = results{2}(:,:,:);
marginals(1,:,:,:) = results{3}(:,:,:);

[C X{i}] = max(marginals, [], 1);
X{i} = reshape(((int32(X{i}).*int32(C > 0.1))), origSize)+1; %Threshold for marginals 
X{i} = permute(X{i}, [1 2 3]);
end
end

function [X,nExamples] = load_1File(imDir)
%Loads IBSR V2 nifti files
nExamples = 1; %for README
X = cell(nExamples);
for i = 1:nExamples
    
    if i < 10
        place = strcat('IBSR_0',int2str(i),'/');
        fileHead = strcat('coreg2IBSR_0',int2str(i));
        %fileHead = strcat('IBSR_0',int2str(i));
    else
        place = strcat('IBSR_',int2str(i),'/');
        fileHead = strcat('coreg2IBSR_',int2str(i));
        %fileHead = strcat('IBSR_',int2str(i));
    end
    
    %Image
    %I_t1uncompress = wfu_uncompress_nifti(strcat(imDir,'coReg2/oppss/',fileHead,'_segTRI_ana.nii'));
    %I_t1uncompress = wfu_uncompress_nifti(strcat(imDir,place,fileHead,'_segTRI_ana.nii'));
    I_t1uncompress = wfu_uncompress_nifti(strcat(imDir,'seg1.hdr'));
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X{i} = I_T1+1;  
    
end
end

function [X,nExamples] = load_fast_file(imDir)
%Loads IBSR V2 nifti files
nExamples = 20; %for README
% X = cell(nExamples);
bList = dir(strcat(imDir,'r*.hdr'));
X = cell(1);
for i = 1:nExamples
    
    %Image
    %I_t1uncompress = wfu_uncompress_nifti(strcat(imDir, 'fast2/IBSR_', num2str(i,'%02d'), '_ana_strip_seg.nii'));
    I_t1uncompress = wfu_uncompress_nifti(strcat(imDir, bList(i).name));
    %I_t1uncompress = wfu_uncompress_nifti(strcat(imDir, 'im', num2str(i,'%02d'), '_seg.nii'));
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    X{i} = I_T1+1;  
%     X{i} = permute(X{i}, [1 3 2]);
%     X{i} = flipdim(X{i},2);
end
for i=1:nExamples
         %1 is background, 2 is CSF, 3 is GM, 4 is WM
         max(max(max(X{i})))
    X{i} = ((X{i}==0) + ((X{i}>120)*2) + ((X{i}>130)*1) + ((X{i}>254)*1)); 
    %X{i} = int32(X{i});
end
end

function [X,y,nExamples] = load_brainweb()
% Currently Only Works for 1 brain. 
if exist('/acmi/chris13/data/brainweb_MS/t1_ai_msles2_3mm_pn3_rf20.rawb', 'file')
  fid = fopen('/acmi/chris13/data/brainweb_MS/t1_ai_msles2_3mm_pn3_rf20.rawb');
  I = fread(fid,inf,'uint8=>int32');
  I = reshape(I, 181, 217, 60);
end

if exist('/acmi/chris13/data/brainweb_MS/phantom_1.0mm_msles2_crisp.rawb', 'file')
  fid = fopen('/acmi/chris13/data/brainweb_MS/phantom_1.0mm_msles2_crisp.rawb');
  y = fread(fid,inf,'uint8=>int32');
  y = reshape(y, 181, 217, 181);
end
I = flipdim(I,2);
y = flipdim(y,2);
y(y>3)=0;
% if exist('/home/cmagnan1/phantom_1.0mm_normal_bck.rawb', 'file')
%   fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_bck.rawb');
%   bck = fread(fid,inf,'uint8=>int32');
%   bck = reshape(bck, 181, 217, 181);
% end
% 
% if exist('/home/cmagnan1/phantom_1.0mm_normal_csf.rawb', 'file')
%   fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_csf.rawb');
%   csf = fread(fid,inf,'uint8=>int32');
%   csf = reshape(csf, 181, 217, 181);
% end
% 
% if exist('/home/cmagnan1/phantom_1.0mm_normal_gry.rawb', 'file')
%   fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_gry.rawb');
%   gry = fread(fid,inf,'uint8=>int32');
%   gry = reshape(gry, 181, 217, 181);
% end
% 
% if exist('/home/cmagnan1/phantom_1.0mm_normal_wht.rawb', 'file')
%   fid = fopen('/home/cmagnan1/phantom_1.0mm_normal_wht.rawb');
%   wht = fread(fid,inf,'uint8=>int32');
%   wht = reshape(wht, 181, 217, 181);
% end


%Create label data
% bck = (bck>=30) * 1;
% csf = (csf>=30) * 2;
% wht = (wht>=30) * 3;
% gry = (gry>=30) * 4;
% y = bck + wht;
% y(y==3) = 2;
% y = y + gry;
% y(y>3) = 3;
% y = y + csf;
% y(y>4) = 4;

X = I;
nExamples = 1;
end

function [X,y,nExamples] = load_ibsr_segs(imDir,segDir)
offsets = [0 1 2 1 1 0 0 0 0 3 3 -4 3 2 6 8 1 0 0 2];
%Load data from ibsr database
bList = dir(strcat(imDir,'*.hdr'));
rawImages = cell(length(bList));
nExamples = length(bList)-1; %CAUSE WHAT THE HELL -4

bList = dir(strcat(segDir,'*.hdr'));
segs = cell(length(bList));
ends = zeros(nExamples,1);

for f = 1:length(bList)
    fid = fopen(strcat(segDir,'r',bList(f).name));
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
end
