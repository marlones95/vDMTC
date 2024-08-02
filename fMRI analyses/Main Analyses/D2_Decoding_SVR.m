%% Decoding script: SVR
% This Script specifies what features should be decoded and then executes the Decoding analysis.
function D2_Decoding_SVR(beta_dir, output_dir, labelnames, dec_type, rad, subj_posterior)


% Create output directory if it does not exist yet
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% =========================================================================
%               Configure input for The Decoding Toolbox
% =========================================================================

clear cfg
cfg = decoding_defaults;
cfg.software = 'SPM12';

% Specify where the results should be saved
cfg.results.overwrite     = 1;
cfg.results.dir           = output_dir;

% DATA SCALING
cfg.scale.method          = 'z';    % standardize the data
cfg.scale.estimation      = 'all';  % only training, only test data, or all

% SEARCHLIGHT SPECIFICATIONS
cfg.analysis              = dec_type;   % = searchlight
cfg.searchlight.unit      = 'voxels';   % Searchlight radius in voxels (not mm) 
cfg.searchlight.radius    = rad;        % 5 voxel radius
cfg.searchlight.spherical = 0;          % only useful for mm
% The amount of information you want to have printed on the screen
% 0: no, 1: normal, output, 2: all):
cfg.verbose               = 0;  

% Method and model parameters 
cfg.decoding.method = 'regression'; % Define SVR as a method

% OUTPUTS SPECIFICATION
cfg.results.output = {'corr', 'zcorr'}; % zcorr is used for the second level

% DISPLAY:
cfg.plot_selected_voxels  = 0; % 0: no plotting, 1: every step, 2: every second step ...
% 100: every hundredth step ...

cfg.plot_design           = 1; % plots the design every time

if dec_type == 'roi'
    cfg.files.mask = masks;
else
    cfg.files.mask = [beta_dir filesep 'mask.nii']; % The mask created with the GLM
end

% =========================================================================
%                       Specify and run the decoding
% =========================================================================

%% Extract subjective frequency differences (SFDs) for all conditions
% First, get the true frequencies for all conditions
stim_set = [16 16 20 20 24 24 28 28;...
            12 20 16 24 20 28 24 32];
% Then, get the subjective log frequencies of the first stimuli (considering
% Fechner's law and the time-order effect)
subj_log_f1            = sort(repmat(unique(subj_posterior.muX(1,:)), 1,2));
% The second frequencies is just the log of the second frequency
% (considering Fechner's law)
subj_log_f2            = log(stim_set(2,:));

% Then calculate the SFDs for all 16 conditions by simple subtraction
subj_log_diffs = zeros(1,16);
for i = 1:8
    subj_log_diffs(i*2-1:i*2)     = [subj_log_f2(i) - subj_log_f1(i), subj_log_f1(i) - subj_log_f2(i)];
end

% Use these SFDs as labels
labels = [subj_log_diffs'];

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

% Extract all information for the cfg.files structure
cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_dir); 


% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg);

% Display the design in the command window
display_design(cfg);


%% Run decoding
results = decoding(cfg);
