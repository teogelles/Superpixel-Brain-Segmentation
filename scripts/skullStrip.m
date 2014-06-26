% Program for stripping the skull off of the ADNI images using
% SPM's tissue segmentations

function skullStrip()
    
    filehead = '/acmi/fmri/';    
    types = {'AD','MCI','CN'};
    totalPatients = [92, 203 ,102];
    fprintf('Stripping skulls ...\n');
    
    for i = 1:3
        type = types{i};
        patientsOfType = totalPatients(i);
        fprintf('Working on %s',type);
        for pnum = 1:2 %patientsOfType
            if mod(pnum,10) == 0
                fprintf('.');
            end
            
            patientfile = strcat('patient',num2str(pnum),'.nii');
            % c1filename = strcat(filehead, type, '_T1/', 'c1', patientfile);
            % c2filename = strcat(filehead, type, '_T1/', 'c2', patientfile);
            % % Since SPM is bad at designating CSF, we decied not to include
            % % since it would also involve including more skull. However, if
            % % one desires otherwise, the code is here
            % c3filename = strcat(filehead, type, '_T1/', 'c3', patientfile);    
            originalfilename = strcat(filehead, type, '_T1/', ...
                                      patientfile);
            chrisfilename = strcat(['/acmi/chris13/results/' ...
                                'ADNIresults/'],type,num2str(pnum),'_again.hdr');
            
            % c1 = load_nifti(c1filename);
            % c2 = load_nifti(c2filename);
            % c3 = load_nifti(c3filename);
            original = load_nifti(originalfilename);
            chris = load_nifti(chrisfilename);
            
            % brain = (c1 ~= 0) + (c2 ~= 0) + (c3 ~= 0);
            % brain = (brain ~= 0);
            brain = (chris ~= 1);
            stripped = original.*brain;
            
            strippedNii = make_nii(stripped);
            % save_nii(strippedNii,strcat('/scratch/tgelles1/summer2014/', ...
            %                             'ADNI_stripped/', type, ...
            %                             num2str(pnum),'.nii'));
            save_nii(strippedNii,strcat('/scratch/tgelles1/summer2014/', ...
                                        'chris_stripped/c', type, ...
                                        num2str(pnum),'.nii'));
            
        end
        fprintf('\n');
    end
    fprintf('Done\n');
end

function ret = load_nifti(filename)
    I_t1uncompress = wfu_uncompress_nifti(filename);
    I_uncompt1 = spm_vol(I_t1uncompress);
    I_T1 = spm_read_vols(I_uncompt1);
    ret = I_T1;
end