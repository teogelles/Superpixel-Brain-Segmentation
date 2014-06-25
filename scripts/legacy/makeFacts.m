%function makeFacts(results, originals, tissues, oName, maxStates)
function makeFacts(fName, typeInd, startI, endI)
%results is a cell array of segmentations
%original is a cell array of the MRI images
%tissues is a cell array of the tissue segmentations(bg=1,csf=2,gm=3.wm=4)
%oName is the file name
    
    oName = fName;
    outF = fopen(fName, 'w');

    [results, originals, tissues] = load_the_things(typeInd, startI, endI);
    nExamples = length(results);
    %maxStates = max(results{startI});
    maxStates = 116;
    for i=startI:endI
        fprintf('Creating features for patient %d\n',i);
        fprintf(outF, '//Patient %d \n', i);
        brain = originals{i};
        segs = results{i};
        tiss = tissues{i};
        has_regions(segs, outF, i, maxStates);
        adj_regions(segs, outF, i, maxStates);
        sizes(segs, outF, i, maxStates);
        centriods(segs, outF, i, maxStates);
        spreads(segs, outF, i, maxStates);
        averages(brain, segs, tiss, outF, i, maxStates);
        var_entrop(brain, segs, outF, i, maxStates);
    end
    close('all');
end

function var_entrop(brain, segs, outF, pid, maxStates)
    for i = 1:maxStates
        variance = var(double(brain(segs==i)));
        fprintf(outF, 'variance(patientid%d, regionid%d, %f).\n', pid, i, variance);
        totalB = sum(sum(sum(segs==i)));
        entropSpread = (hist(double(brain(segs==i)),255)./totalB).*log(hist(double(brain(segs==i)),255)./totalB);
        entrop = -1 * sum(entropSpread(~isnan(entropSpread))); 
        fprintf(outF, 'entropy(patientid%d, regionid%d, %f).\n', pid, i, entrop);
    end
end

function averages(brain, segs, tissues, outF, pid, maxStates)
    size(brain)
    size(tissues)
    size(segs)
    for i = 1:maxStates
        totalAvg = mean(brain(segs==i));
        t = brain .* int32(tissues == 4);
        wmAvg = mean(t(segs==i & tissues==4));
        t = brain .* int32(tissues == 3);
        gmAvg = mean(t(segs==i & tissues==3));
        t = brain .* int32(tissues == 2);
        csfAvg = mean(t(segs==i & tissues==2));
        fprintf(outF, 'avgintensity(patientid%d, regionid%d, %f).\n',pid, i, totalAvg);
        fprintf(outF, 'wmavg(patientid%d, regionid%d, %f).\n',pid, i, wmAvg);
        fprintf(outF, 'gmavg(patientid%d, regionid%d, %f).\n',pid, i, gmAvg);
        fprintf(outF, 'csfavg(patientid%d, regionid%d, %f).\n',pid, i, csfAvg);
    end
end

function spreads(segs, outF, pid, maxStates)
    for i=1:maxStates
        idx = find(segs == i);
        [rows cols pags] = ind2sub(size(segs),idx);
        spreadX = (max(rows)-min(rows))/min(rows);
        spreadY = (max(cols)-min(cols))/min(cols);
        spreadZ = (max(pags)-min(pags))/min(pags);
        fprintf(outF, 'spreadx(patientid%d, regionid%d, %d).\n', pid, i, int32(spreadX));
        fprintf(outF, 'spready(patientid%d, regionid%d, %d).\n', pid, i, int32(spreadY));
        fprintf(outF, 'spreadz(patientid%d, regionid%d, %d).\n', pid, i, int32(spreadZ));
    end
end

function has_regions(segs, outF, pid, maxStates)
    maxStates
    for i=1:maxStates
        if(any(any(any(segs == i))))
            fprintf(outF,'has_region(patientid%d, regionid%d).\n',pid,i);
        end
    end
end

function adj_regions(segs, outF, pid, maxStates)
%FIND SOME WAY WITHOUT LOOPING OVER EACH
%STRFIND??? COVAR???
end

function sizes(segs, outF, pid, maxStates)
    for i=1:maxStates
        size = sum(sum(sum(segs==i)));
        if (size>0)
            fprintf(outF, 'size(patientid%d,regionid%d,%d).\n', pid,i,size);
        end
    end
end

function centriods(segs, outF, pid, maxStates)

    for i=1:maxStates
        if (any(any(any(segs==i))))
            M = (segs == i);
            [rows cols slices] = size(M);
            
            y = 1:rows;
            x = 1:cols;
            z = 1:slices;
            
            [X, Y, Z] = meshgrid(x,y,z);
            
            cY = mean(Y(M));
            cX = mean(X(M));
            cZ = mean(Z(M));
            
            fprintf(outF, 'centroidx(patientid%d, regionid%d, %f).\n',pid,i,cX);
            fprintf(outF, 'centroidy(patientid%d, regionid%d, %f).\n',pid,i,cY);
            fprintf(outF, 'centroidz(patientid%d, regionid%d, %f).\n',pid,i,cZ);
        end
    end
end

function [anaSegs, originals, tissues] = load_the_things(typeInd, startI, endI)

    types = {'MCI' 'AD' 'CN'};
    nums = [0 0 0];

    %load the originals
    originals = cell(1,1);
    for i = typeInd:typeInd
        origDir = strcat('/acmi/fmri/',types{i}, '_T1/');
        bList = dir(strcat(origDir,'patient*.nii')); %constructs list of wanted files
        nIms = length(bList);
        nums(i) = nIms;
        nums(i) = endI; nIms = endI; %AHHH TESTING CODE
        tempIms = cell(nIms,1);
        for j=startI:nIms
            strcat(origDir,'patient',num2str(j),'.nii')
            I_t1uncompress = wfu_uncompress_nifti(strcat(origDir,'patient',num2str(j),'.nii'));
            I_uncompt1 = spm_vol(I_t1uncompress);
            I_T1 = spm_read_vols(I_uncompt1);
            tempIms{j} = int32(I_T1);
        end
        originals{1} = tempIms;
    end
    types = {'Mci' 'Ad' 'Cn'};
    nums
    %load the results
    anaSegs = cell(1,1);

    for i = typeInd:typeInd
        anaDir = strcat('/acmi/fmri/altAtlas/');
        bList = dir(strcat(anaDir, strcat(types{i},'*CO*.nii'))); %constructs list of wanted files
        nIms = nums(i);
        tempIms = cell(nIms,1);
        for j=startI:nIms
            I_t1uncompress = wfu_uncompress_nifti(strcat(anaDir,'CR',types{i},'At',num2str(j),'aal_MNI_V4.nii'));
            I_uncompt1 = spm_vol(I_t1uncompress);
            I_T1 = spm_read_vols(I_uncompt1);
            tempIms{j} = int32(I_T1);
        end
        anaSegs{1} = tempIms;
    end
    types = {'MCI' 'AD' 'CN'};
    %load the tissue segmentations
    tissues = cell(3,1);
    for i = typeInd:typeInd
        tissDir = strcat('/acmi/chris13/results/ADNIresults/');
        nIms = nums(i);
        tempIms = cell(nIms,1);
        for j=startI:nIms
            %bList = dir(strcat(tissDir,types{i},num2str(j),'*'));
            %indx = length(bList);
            %I_t1uncompress = wfu_uncompress_nifti(strcat(tissDir,bList(indx).name));
            I_t1uncompress = wfu_uncompress_nifti(strcat(tissDir,types{i},num2str(j),'_again.hdr'));
            I_uncompt1 = spm_vol(I_t1uncompress);
            I_T1 = spm_read_vols(I_uncompt1);
            tempIms{j} = int32(I_T1);
        end
        tissues{1} = tempIms;
    end

    tissues = flatten(tissues);
    originals = flatten(originals);
    anaSegs = flatten(anaSegs);


end


function y = flatten(x)
% y = flatten(x)
%
% Takes a cell array x that contains nested cell arrays and
% flattens the contents into a single cell array.
% E.g. flatten({1, {2, 3, {4}, 5}}) returns {1, 2, 3, 4, 5}

% Copyright (C) 2007 Ron J. Weiss
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    if ~iscell(x)
        error('flatten only works on cell arrays.');
    end

    y = inner_flatten(x);
end

function y = inner_flatten(x)
    if ~iscell(x)
        y = {x};
    else
        y = {};
        for n = 1:length(x)
            tmp = inner_flatten(x{n});
            y = {y{:} tmp{:}};
        end
    end
end
