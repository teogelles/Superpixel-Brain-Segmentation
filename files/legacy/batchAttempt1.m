function batchAttempt1()
% List of open inputs
nrun = 203; % enter the number of runs here
jobfile = {'/acmi/chris13/scripts/batchCOREG_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
    theOtherJob = the_other_job(crun);
    spm('defaults', 'FMRI');
    spm_jobman('serial', theOtherJob, '', inputs{:});
end

% List of open inputs
%nrun = 2; % enter the number of runs here
% jobfile = {'/acmi/chris13/scripts/batchAttempt1_job.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(0, nrun);
% for crun = 1:nrun
%     theJob = the_job(crun);
%     spm('defaults', 'FMRI');
%     spm_jobman('serial', theJob, '', inputs{:});
% end

end

function matlabbatch = the_other_job(b)
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.coreg.write.ref = {strcat('/acmi/fmri/MCI_T1/patient',num2str(b),'.nii')};
matlabbatch{1}.spm.spatial.coreg.write.source = {strcat('/acmi/fmri/altAtlas/MciAt',num2str(b),'aal_MNI_V4.nii')};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix ='CR';


end


function matlabbatch = the_job(b)
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {strcat('/acmi/fmri/CN_T1/patient',num2str(b),'_seg_inv_sn.mat')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {strcat('/acmi/fmri/altAtlas/CNCO',num2str(b),'aal_MNI_V4.nii,1')};
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = strcat('CNAT',num2str(b));
end
