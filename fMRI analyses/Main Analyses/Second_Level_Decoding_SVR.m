%% Second level analysis performing a one-sample t-test, to statistically test for local brain activation patterns
clc; clear;

% =========================================================================
%                   Specify parameters and analyses
% =========================================================================

data_dir        = 'D:\vDMTC\data';
GLM_name        = 'Decoding_SVR';
radius          = 5;
analysis_name   = ['Searchlight\' 'SVM_r' num2str(radius) 'vox\192s\' GLM_name];


% sub-025 excluded because of problems with projector
SJs     = {'sub-001','sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008',...
    'sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024',...
    'sub-026','sub-027','sub-028','sub-029','sub-030','sub-031','sub-032','sub-033','sub-034','sub-035',...
    'sub-036','sub-037'};

group_mask      = '36Subjects_SVR'; % mean MNI-space mask for the second level Analysis from the first level masks

% Output directory:
out_dir         = fullfile(data_dir, 'RFX\TTest', [analysis_name]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
cd(out_dir)

% =========================================================================
%                   Select the decoding files
% =========================================================================

fs_all=[];
for sj = 1:length(SJs)

    % select the smoothed, normalized images
    filt    = ['^s3w'];
    % Source directory for the 'mean'-images
    sj_dir =  fullfile(data_dir, SJs{sj}, 'dec_SVR');
    % select the files
    f = spm_select('FPList', sj_dir, filt);
    numVols = size(f,1);
    % create SPM style file list for model specification
    fs_all = [fs_all; cellstr(strcat(f, ',1'))];
end

% =========================================================================
%                       Setting up the matlabbatch
% =========================================================================

% create SPM style file list for model specification
matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(fs_all);
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = ...
    cellstr(fullfile(data_dir, 'MNIMask', ['GroupMask_' group_mask '.nii']));
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run',matlabbatch);
clear matlabbatch


% =========================================================================
%               Model estimation and contrast specification
% =========================================================================

%% Model estimation
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[pwd filesep 'SPM.mat']}; % pwd = current location
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run',matlabbatch);
clear matlabbatch;

%% Contrast specification
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[pwd filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = ['pos'];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name       = ['neg'];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec     = -1;
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep    = 'none';
matlabbatch{1}.spm.stats.con.delete                     = 0;

spm_jobman('run',matlabbatch);
clear matlabbatch;