%% Script setting up the the conjunction analysis

% =========================================================================
%                           Run the t-test
% =========================================================================
nrun = 1; % enter the number of runs here
jobfile = {'D:\vDMTC\fMRI Analysis\Conjunction\t_test_job_SVM.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

% =========================================================================
%                   Estimating the GLM Parameters
% =========================================================================
nrun = 1; % enter the number of runs here
jobfile = {'D:\vDMTC\fMRI Analysis\Conjunction\Estimate_job_SVM.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch;


% =========================================================================
%                   Defining and running contrasts
% =========================================================================
n1 = 36; % sample size of my data
n2 = 30; % sample size of Wu (2021)
matlabbatch{1}.spm.stats.con.spmmat(1)                  = ...
    {['D:\vDMTC\data\RFX\Conjunction\SVM' filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = ['vDMTC'];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = [1 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name       = ['SFC-II'];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec     = [0 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep    = 'none';
matlabbatch{1}.spm.stats.con.delete                     = 0;

spm_jobman('run',matlabbatch);
clear matlabbatch;