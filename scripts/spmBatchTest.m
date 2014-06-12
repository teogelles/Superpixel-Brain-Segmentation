% List of open inputs
nrun = 5; % enter the number of runs here
jobfile = {'/acmi/chris13/scripts/spmBatchTest_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
    spm('defaults', 'FMRI');
    spm_job(fName)
    spm_jobman('serial', jobs, '', inputs{:});
end

function spm_job(fName)
matlabbatch{1}.spm.spatial.realign.estwrite.data = {
                                                    {
                                                    '/acmi/fmri/CN_T1/stupidExamplempatient35.nii,1'
                                                    '/acmi/fmri/CN_T1/stupidExamplempatient36.nii,1'
                                                    '/acmi/fmri/CN_T1/stupidExamplempatient37.nii,1'
                                                    '/acmi/fmri/CN_T1/stupidExamplempatient38.nii,1'
                                                    }
                                                    }';
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
matlabbatch{2}.spm.spatial.coreg.estwrite.ref = {'/acmi/fmri/CN_T1/patient99.nii,1'};
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1) = fName;
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'image';
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).sname = 'Realign: Estimate & Reslice: Mean Image';
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','rmean');
matlabbatch{2}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
end