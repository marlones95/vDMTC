%% MVPA support vector machine classification batch script

% This Batch Script first specifies what features should be decoded and
% executes it.

clc; clear;

% =========================================================================
%                   Specify parameters and analyses
% =========================================================================

% Add vector of subject numbers (e.g. 1,2,3,4,5) if any subject has been
% already analysed to exclude that subjects from the analysis
excludeSJ = [];

% Specify the current prefix (here ra because the slice-time corrected and
% realigned images are used)
currPrefix = ['ra'];


first_level_glm = 1; % First level GLM
dec_anal        = 1; % Support vector regression analysis
coregister      = 1; % Coregister functional images to T1
segmentation    = 1; % Segmentation of T1
normalization   = 1; % Normalization
smoothing       = 1; % Smoothing

% Add necessary paths
addpath('C:\Users\nnu16\Documents\MATLAB\spm12'); % SPM
addpath('C:\Users\nnu16\Documents\MATLAB\decoding_toolbox'); % The Decoding Toolbox
addpath('D:\vDMTC\fMRI Analysis\Preprocessing'); % Path containing pre-processing scripts
addpath('D:\vDMTC\fMRI Analysis\Decoding'); % Path containing decoding scripts

% Specify data paths
targetF = 'D:\vDMTC\data';    % target directory for all nifti files

% =========================================================================
%           Step 1: Data organization and analysis configuration
% =========================================================================

%%% Data organization %%%
A2_data_org
% A2_data_org does the following:
% - create a cell array Subjects with Subject names
% - unzip all nii.gz files
% - create a cell array containing filenames for all runs per subject
% - compares info from json-file to niftimetadata and gives warning if they do not match
% - gets the TR and the number of slices from the json file


%%% Analysis configuration %%%
dec_type = "searchlight"; % specifies to perform searchlight decoding

% Choose what you want to decode: '_left_vs_right' for motor response, '_rule' for rule decoding
version = '_left_vs_right'; 
if strcmp('_left_vs_right',version)
    condnames = {'decision_left', 'decision_right'}; 
elseif strcmp('_rule',version) || strcmp('_rule_sub',version)
    condnames = {'FAS', 'SAF'}; % First against second, Second against first
end

duration = 0; % epoch duration, 0 because it is a single event
refslice=slice_order(round(length(slice_order)/2)); % reference slice for slice-time correction
hpf      = 192; % High-pass filter cut-off;
hm = 1;   % include head motion parameters from realignment
CSFWM_params = 1;  % include the first five principal components explaining variance
% in the white matter and cerebrospinal fluid signals, respectively


% =========================================================================
%                           Step 2: First level GLM
% =========================================================================

if first_level_glm
    outputfolder = append('\decoding_glm',version);
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            display(['Step 2, 1st level glm: ' Subjects{sj} ])
            subj_dir = fullfile(targetF, Subjects{sj});
            % D1_glm_1stLevel_left_vs_right: GLM for motor response
            % D1_glm_1stLevel_rule: GLM for rule
            if strcmp('_left_vs_right',version)
                D1_glm_1stLevel_left_vs_right(Subjects, sj, outputfolder, TR, hpf, runs, duration, hm, CSFWM_params);
            elseif strcmp('_rule',version)
                D1_glm_1stLevel_rule(Subjects, sj, outputfolder, TR, hpf, runs, duration, hm, CSFWM_params);
            end
        end
    end
end

% =========================================================================
%                           Step 3: Decoding
% =========================================================================

if dec_anal
    labelnames = condnames;
    rad = 5;
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            display(['Step 3, Decoding: ' Subjects{sj} ])
            beta_dir = fullfile(targetF, Subjects{sj}, append('decoding_glm',version));
            output_dir = fullfile(targetF, Subjects{sj}, append('dec',version));
            % Performs decoding for the respective version
            D2class_Decoding(beta_dir, output_dir, labelnames, dec_type, rad);
        end
    end
end

% =========================================================================
%                       Step 4: Coregistration
% =========================================================================

if coregister
    Co_er = 0; % default setting: only estimate (no reslice), if 1, then estimate & reslice
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            if exist([targetF filesep Subjects{sj} filesep 'func'])
                funcPath    = [targetF filesep Subjects{sj}];
                func_dir    = fullfile(funcPath, 'func');
                struct_dir  = fullfile(funcPath, 'anat');
                sj_dir      = [targetF filesep Subjects{sj}];
                data_dir    = [sj_dir filesep append('dec',version)];
                f3          = spm_select('List', data_dir, '^res_accuracy_minus_chance.nii'); % select accuracy-chance file
                numVols     = size(f3,1);
                Images      = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
                if Co_er ~= 1
                    display(['Step 3, coregistration (estimate): ' Subjects{sj}])
                    % A5a_coregister_est performs coregistration using the 
                    % spm_jobman with the pre-defined parameters. It 
                    % coregisters the structural image to the mean of the 
                    % functional images. Further, the decoding image 
                    % (res_zcorr.nii) is also coregistered.
                    A5a_coregister_est(currPrefix, func_dir, struct_dir, sj, '^sub.*\.nii', Images);
                else
                    display(['Step 5b, coregistration (estimate & reslice): ' Subjects{sj}])
                    A5b_coregister_est_re(currPrefix, func_dir, struct_dir, Subjects{sj}, '^sub.*\.nii', Images);
                end
            else
                display('###########################################################')
                display(['############### ' Subjects{sj} 's functional data do not exist ###########'])
            end
        end
    end
end

% =========================================================================
%                           Step 5: Segmentation
% =========================================================================

% Segmentation --> prefix y_ (to structural image)
if segmentation
    warning off
    SPM_path  = 'C:\Users\nnu16\Documents\MATLAB\spm12';
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        elseif exist([targetF filesep Subjects{sj} filesep ana])==7
            display(['Step 8, segmentation: ' Subjects{sj}])
            struct_dir = fullfile(targetF, Subjects{sj}, 'anat');
            % A6_segmentation performs segmentation of the structural image
            % using the spm_jobman with pre-defined parameters. 
            A6_segmentation(struct_dir, Subjects{sj}, SPM_path, '^sub.*\.nii');
        else
            display('###########################################################')
            display(['############### ' Subjects{sj} ', ' ana ' does not exist ###########'])
        end
    end
end


% =========================================================================
%                           Step 6: Normalization
% =========================================================================

% normalization --> prefix: w
if normalization
    vox_size = [2 2 2]; % Normalize to 2x2x2 voxel size
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            if exist([targetF filesep Subjects{sj} filesep append('dec',version)])
                display(['Step 6, normalization: ' Subjects{sj}])
                funcPath = [targetF filesep Subjects{sj}];
                struct_dir = fullfile(funcPath, 'anat');
                data_dir = fullfile(funcPath, append('dec',version));
                % A7_normalization_SVR performs normalization of functional
                % image to MNI space using the spm_jobman with pre-defined parameters.
                A7_normalization_SVR(data_dir, struct_dir, vox_size);
            else
                display('###########################################################')
                display(['############### ' Subjects{sj} ', ' runs{sj, r} ' does not exist ###########'])
            end
        end
    end
end

% =========================================================================
%                           Step 7: Smoothing
% =========================================================================

% smoothing --> prefix: s
if smoothing
    kernel_size = [3 3 3]; % smoothing kernel of 3mm

    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            display(['Step 7, Smoothing: ' Subjects{sj}])
            sj_dir = [targetF filesep Subjects{sj}];
            data_dir = [sj_dir filesep append('dec',version)];
            % A8_smoothing_run performs smoothing of functional image with
            % a pre-defined kernel using the spm_jobman with pre-defined parameters.
            A8_smoothing_run(data_dir, Subjects{sj}, ['^wres_accuracy_minus_chance.nii'], kernel_size);
        end
    end
end