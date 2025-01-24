%% Script setting up the t-test for the conjunction
matlabbatch{1}.spm.stats.factorial_design.dir = {'D:\vDMTC\data\RFX\Conjunction\SVM'};

% Loading in all 1st level prediction accuracy maps from both studies
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {
                                                           'D:\vDMTC\data\sub-002\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-003\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-004\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-005\dec\s3wres_accuracy_minus_chance.nii,1'  
                                                           'D:\vDMTC\data\sub-007\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-008\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-009\dec\s3wres_accuracy_minus_chance.nii,1'  
                                                           'D:\vDMTC\data\sub-010\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-012\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-016\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-017\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-018\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-022\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-023\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-026\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-027\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-029\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-030\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-031\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-032\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-033\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-034\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-036\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\sub-037\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           };
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = {
                                                           'D:\vDMTC\data\YH\sub-001\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-002\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-003\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-004\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-005\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-006\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-007\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-008\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-009\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-010\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-011\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-012\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-013\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-014\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-015\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-016\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-017\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-018\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-019\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-020\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-021\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-022\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-023\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-024\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-025\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-026\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           'D:\vDMTC\data\YH\sub-027\dec\s3wres_accuracy_minus_chance.nii,1'
                                                           };
% Further specifications
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = ...
    {'D:\vDMTC\fMRI Analysis\full_brain_mask.nii'}; % Full brain mask
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
