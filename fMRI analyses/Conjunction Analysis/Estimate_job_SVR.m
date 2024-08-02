%% Script setting up the parameter estimation of the conjunction
matlabbatch{1}.spm.stats.fmri_est.spmmat = {['D:\vDMTC\data\RFX\Conjunction\SVR\SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
