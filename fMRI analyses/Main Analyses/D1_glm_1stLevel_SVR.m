%% First Level GLM for Decoding Subjective Frequency Differences

% Setting up the GLM for the Support vector regression by estimating betas for all 16 conditions per run
function D1_glm_1stLevel_SVR(Subjects, sj, outputfolder_1st,...
    tr, hpf, runs, condnames, duration, hm, CSFWM_params)


% =========================================================================
%                           Initialization
% =========================================================================

% load the log_file
log_dir = 'D:\vDMTC\Logs';
current_sj = Subjects{sj}(end-1:end);
logfile = spm_select('FPList', log_dir, ['^fMRI_', current_sj, '.*\.mat$']);
load(logfile);


% set SPM defaults
spm('defaults','fmri')
spm_jobman('initcfg');

% data directory of the current_sj
data_dir = ['D:\vDMTC\data\sub-0', current_sj];

% target directory that will contain the created job.mat file and SPM file
tgt_dir = [data_dir, outputfolder_1st];

% directory with movement parameters and functional runs
func_dir = [data_dir, '\func'];

if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end
cd(tgt_dir)

% =========================================================================
%                           Extract timings
% =========================================================================

% Get timing of delayperiod for all trials
delayperiod = squeeze(mylog.timing(5,:,:))';

%% Extract timings of all conditions
% Initialize structures
idx = struct();         % Indexing trials
decision = struct();    % Timing of all conditions

for i = 1:length(condnames)
    currcond        = condnames{i};               % condition name
    f1              = str2double(currcond(1:2));  % f1 frequency
    f2              = str2double(currcond(4:5));  % f2 frequency
    condition_type  = currcond(7:end);            % defines the correct decision (high or low)
    
    if (strcmp(condition_type, 'low') && f1 < f2) || ...
       (strcmp(condition_type, 'high') && f1 > f2)
        foc = 1; % rule = compare f1 vs f2
    elseif (strcmp(condition_type, 'low') && f1 > f2) || ...
           (strcmp(condition_type, 'high') && f1 < f2)
        foc = 2; % rule = compare f2 vs f1
    end

    % Create dynamic field names
    field_name = matlab.lang.makeValidName(currcond); % Ensure valid field name
    
    % Compute logical indexing
    idx.(field_name) = (mylog.focus == foc) & ...
                       (mylog.flutter(:,:,1) == f1) & ...
                       (mylog.flutter(:,:,2) == f2);
    % mylog.focus indicates which stimulus is the baseline (2 = f1 is baseline; 1 = f2 is Baseline) --> 
    % if focus == 1 & mylog.flutter(1) > mylog.flutter(2) or if focus == 2 & mylog.flutter(1) < mylog.flutter(2)
    % --> decision high and vice versa

    % Get timings for the current condition
    decision.(field_name) = zeros(6, 64);
    decision.(field_name)(idx.(field_name)) = delayperiod(idx.(field_name));

end

% =========================================================================
%                       Setting up the matlabbatch
% =========================================================================

% event-related analysis
matlabbatch{1, 1}.spm.stats.fmri_spec.dir = cellstr(tgt_dir); % Output Directory
matlabbatch{1, 1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1, 1}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{1, 1}.spm.stats.fmri_spec.timing.fmri_t = 16; % Temporal resolution
matlabbatch{1, 1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % reference slice

filter = ['^rasub-0', current_sj, '_task-task_run-'];
f = spm_select('List', func_dir, filter);

for r = 1:size(runs,2)
    V=spm_vol([func_dir filesep f(r,:)]);
    files={};
    for i=1:(size(V,1))
        files{i} = [func_dir filesep strtrim(f(r,:)) ',' int2str(i)];
    end    


    % Allocatin the conditions
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).scans = cellstr(files');

    matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf; % high pass filter

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(1).name      = '16_12_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(1).onset     = nonzeros(decision.x16_12_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(1).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(2).name      = '16_12_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(2).onset     = nonzeros(decision.x16_12_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(2).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(3).name      = '16_20_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(3).onset     = nonzeros(decision.x16_20_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(3).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(4).name      = '16_20_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(4).onset     = nonzeros(decision.x16_20_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(4).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(5).name      = '20_16_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(5).onset     = nonzeros(decision.x20_16_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(5).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(6).name      = '20_16_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(6).onset     = nonzeros(decision.x20_16_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(6).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(7).name      = '24_20_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(7).onset     = nonzeros(decision.x24_20_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(7).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(8).name      = '24_20_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(8).onset     = nonzeros(decision.x24_20_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(8).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(9).name      = '20_24_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(9).onset     = nonzeros(decision.x20_24_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(9).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(10).name      = '20_24_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(10).onset     = nonzeros(decision.x20_24_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(10).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(11).name      = '24_28_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(11).onset     = nonzeros(decision.x24_28_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(11).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(12).name      = '24_28_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(12).onset     = nonzeros(decision.x24_28_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(12).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(13).name      = '28_24_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(13).onset     = nonzeros(decision.x28_24_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(13).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(14).name      = '28_24_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(14).onset     = nonzeros(decision.x28_24_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(14).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(15).name      = '28_32_low';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(15).onset     = nonzeros(decision.x28_32_low(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(15).duration  = duration;

    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(16).name      = '28_32_high';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(16).onset     = nonzeros(decision.x28_32_high(r,:))';
    matlabbatch{1, 1}.spm.stats.fmri_spec.sess(r).cond(16).duration  = duration;

    % Motion Parameters and CSF
    if hm && CSFWM_params == 0
        mf = spm_select('FPList', func_dir, ['rp_asub-0', current_sj, '_task-task_run-', num2str(r)]);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {mf};
    elseif hm && CSFWM_params
        mf = spm_select('FPList', func_dir, ['rp_asub-0', current_sj, '_task-task_run-', num2str(r)]);
        CSF = spm_select('FPList', [data_dir '\PCA_4mm9090'], ['^run0' num2str(r) '_csf.mat']);
        WM = spm_select('FPList', [data_dir '\PCA_4mm9090' ], ['^run0' num2str(r) '_wm.mat']);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {CSF; WM; mf};
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf; % high pass filter

end

matlabbatch{1, 1}.spm.stats.fmri_spec.fact  	        = struct('name', {}, 'levels', {});
% model hrf and first temporal derivative
matlabbatch{1, 1}.spm.stats.fmri_spec.bases.hrf .derivs  = [0 0];
% model interactions
matlabbatch{1, 1}.spm.stats.fmri_spec.volt              = 1;
% global normalization
matlabbatch{1, 1}.spm.stats.fmri_spec.global            = 'None';
% masking threshold
matlabbatch{1, 1}.spm.stats.fmri_spec.mthresh           = 0.8;
% explicit mask
matlabbatch{1, 1}.spm.stats.fmri_spec.mask              = {''};
% autocorrelation modelling (whitening filter)
matlabbatch{1, 1}.spm.stats.fmri_spec.cvi               = 'AR(1)';


% =========================================================================
%                               Run the GLM
% =========================================================================

%% Create the model
fprintf('Creating GLM\n')
display(['Univariate First Level: SJ' current_sj])
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch

%%  Model Estimation:
load(fullfile(tgt_dir, filesep, 'SPM.mat'));
fprintf('Estimating GLM \n');
cd(tgt_dir);
SPM = spm_spm(SPM);

clear SPM;
cd(log_dir);