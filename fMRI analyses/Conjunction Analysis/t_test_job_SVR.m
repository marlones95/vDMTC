%% Script setting up the t-test for the conjunction
matlabbatch{1}.spm.stats.factorial_design.dir = {'D:\vDMTC\data\RFX\Conjunction\SVR'};

% Loading in all 1st level prediction accuracy maps from both studies
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {
                                                           'D:\vDMTC\data\sub-001\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-002\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-003\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-004\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-005\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-006\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-007\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-008\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-009\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-010\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-011\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-012\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-013\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-014\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-015\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-016\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-017\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-018\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-019\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-021\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-022\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-023\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-024\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-026\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-027\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-028\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-029\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-030\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-031\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-032\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-033\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-034\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-035\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-036\dec_SVR\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\data\sub-037\dec_SVR\s3wres_zcorr.nii,1'
                                                           };
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = {
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-001\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-002\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-003\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-004\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-005\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-006\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-007\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-008\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-009\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-010\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-011\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-012\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-013\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-014\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-015\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-016\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-017\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-018\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-019\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-020\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-021\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-022\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-023\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-024\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-025\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-026\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-027\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-028\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-029\dec_SVR_Wu\s3wres_zcorr.nii,1'
                                                           'D:\vDMTC\Yuan Hao_My_Analyses\data\sub-030\dec_SVR_Wu\s3wres_zcorr.nii,1'
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
    {'D:\vDMTC\fMRI Analysis\Conjunction\MNIMask\GroupMask_ConjunctionMask66Subjects.nii,1'}; % Conjunction mask from all subjects
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
