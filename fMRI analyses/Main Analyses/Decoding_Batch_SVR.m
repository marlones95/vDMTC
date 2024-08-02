%% MVPA support vector regression analysis batch script

% This Batch Script specifies preprocessing and decoding analyses and
% executes the respective scripts

clc; close all; clear all;

% =========================================================================
%                   Specify parameters and analyses
% =========================================================================

% Add vector of subject numbers (e.g. 1,2,3,4,5) if any subject has been
% already analysed to exclude that subjects from the analysis
excludeSJ = [];

% Specify the current prefix (e.g. ra if the slice-time corrected and
% realigned images should be used)
currPrefix = [''];

dicom_import    = 1; % Import dicom to nifti
slicetiming     = 1; % Slice-time correction
realignment     = 1; % realigment
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
sourceF = 'D:\vDMTC\rawdata'; % source directory containing raw dicom files
targetF = 'D:\vDMTC\data';    % target directory for all nifti files


% =========================================================================
%       Step 1: Import dicom files to BIDS formatted nifti files
% =========================================================================

if dicom_import
    sess = [0 1]; % First entry: Do wou want an extra instance for session folders?
    % (yes = 1, no = 0)
    % Second entry: How many sessions are there per subject? (eg. 1 or 4) Not runs!
    if ~exist(targetF, 'dir')
        mkdir(targetF);
    end

    cd(sourceF);
    pb = dir('*ccnb_*'); % List of all directories with the relevant raw data
    task = 'task';

    for subj = 1:length(pb)
        if ismember(subj, excludeSJ)
            continue;
        else
            if subj < 10
                sub = ['00' num2str(subj)];
            elseif subj < 100
                sub = ['0' num2str(subj)];
            else
                sub = num2str(subj);
            end
            fSubj = [targetF '/sub-' sub '/']; % Name subject directory
            if ~exist(fSubj, 'dir'), mkdir(fSubj); end
            dcmDir = [pb(subj).folder filesep pb(subj).name]; % Corresponding dicom directory
            % Call dicm2bids function that uses Xiangrui Li's dicm2nii.m
            A1_dicm2bids(dcmDir, targetF, subj, task, sess);
        end
    end
end

% =========================================================================
%           Step 2: Data organization and analysis configuration
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

% Specify condition names (first number refers to f1 frequency, second to
% f2 frequency, low or high refers to the correct decision. E.g. in
% 16_12_low participants had to compare f2 against f1, such that the
% correct decision was lower, while in 16_12_high participants had to
% compare f1 against f2 so the correct decision was higher
condnames = {'16_12_low','16_12_high','16_20_low','16_20_high','20_16_low','20_16_high',...
    '20_24_low','20_24_high','24_20_low','24_20_high','24_28_low','24_28_high',...
    '28_24_low','28_24_high','28_32_low','28_32_high'};
duration = 0; % epoch duration, 0 because it is a single event
refslice=slice_order(round(length(slice_order)/2)); % reference slice for slice-time correction
hpf      = 192; % High-pass filter cut-off;
hm = 1;   % include head motion parameters from realignment
CSFWM_params = 1;  % include the first five principal components explaining variance
% in the white matter and cerebrospinal fluid signals, respectively

% =========================================================================
%                       Step 3: Slice-time correction
% =========================================================================

% slice-timing --> prefix: a
if slicetiming
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            for r = 1:size(runs, 2)
                if exist([targetF filesep Subjects{sj} filesep 'func' filesep runs{sj, r}])
                    display(['Step 3, slice time correction: ' Subjects{sj} ', ' runs{sj, r}])
                    funcPath = [targetF filesep Subjects{sj} filesep 'func'];
                    run_dir = fullfile(funcPath);
                    % A3_slice_time_correction performs slice-time correction
                    % using the spm_jobman with the pre-defined parameters
                    A3_slice_time_correction(Subjects{sj},runs{sj, r}, run_dir,...
                        ['^' currPrefix runs{sj, r}],n_slices,slice_order,refslice,TR);
                else
                    display('###########################################################')
                    display(['############### ' Subjects{sj} ', '...
                        runs{sj, r} ' does not exist ###########'])
                end
            end
        end
    end
    currPrefix=['a' currPrefix]; % add the prefix a
end

% =========================================================================
%                       Step 4: Realignment
% =========================================================================

% Realignment --> prefix: r
if realignment
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            sj_dir = [targetF filesep Subjects{sj}];
            if exist([targetF filesep Subjects{sj} filesep 'func' filesep runs{sj, r}])
                display(['Step 4, realignment: ' Subjects{sj} ', ' runs{sj, r}])
                funcPath = [targetF filesep Subjects{sj} filesep 'func'];
                run_dir = fullfile(funcPath);
                for r = 1:size(runs, 2)
                    run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                end
                % A4_realignment performs Realignment: Estimate & Reslice 
                % using the spm_jobman with the pre-defined parameters
                A4_Realignment([sj_dir filesep 'func'], run_files);
            else
                display('###########################################################')
                display(['############### ' Subjects{sj} ', '...
                    runs{sj, r} ' does not exist ###########'])
            end
        end

    end
    currPrefix=['r' currPrefix]; % add the prefix r
end

% =========================================================================
%                           Step 5: First level GLM
% =========================================================================

if first_level_glm
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            outputfolder = '\decoding_glm_SVR';
            display(['Step 5, 1st level glm: ' Subjects{sj}])
            subj_dir = fullfile(targetF, Subjects{sj});
            % D1_glm_1stLevel_SVR sets up GLM estimating betas for all 16
            % conditions
            D1_glm_1stLevel_SVR(Subjects, sj, outputfolder, TR,...
                hpf, runs, condnames, duration, hm, CSFWM_params);
        end
    end
end

% =========================================================================
%                           Step 6: Decoding
% =========================================================================

if dec_anal
    labelnames = condnames;
    rad = 5; % Searchlight radius
    % Get the posteriors from the Bayesian bevioural model
    load('D:\vDMTC\Behavioural Analysis\Bayesian Model\VBA\posteriors.mat');
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            display(['Step 6, Decoding: ' Subjects{sj} ])
            beta_dir = fullfile(targetF, Subjects{sj}, 'decoding_glm_SVR');
            output_dir = fullfile(targetF, Subjects{sj}, 'dec_SVR');
            % D2regress_Decoding_batch runs the SVR analysis with the SFDs
            % as labels           
            D2_Decoding_SVR(beta_dir, output_dir, labelnames, dec_type, rad, ...
                subj_posterior(sj));
        end
    end
end

% =========================================================================
%                       Step 7: Coregistration
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
                data_dir    = [funcPath filesep 'dec_SVR'];
                f3          = spm_select('List', data_dir, '^res_zcorr.nii'); % select zcorr file created by the SVR
                numVols     = size(f3,1);
                Images      = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
                if Co_er ~= 1
                    display(['Step 7, coregistration (estimate): ' Subjects{sj}])
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
%                           Step 8: Segmentation
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
%                           Step 9: Normalization
% =========================================================================

% normalization --> prefix: w
if normalization
    vox_size = [2 2 2]; % Normalize to 2x2x2 voxel size
    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            if exist([targetF filesep Subjects{sj} '\dec_SVR'])
                display(['Step 9, normalization: ' Subjects{sj}])
                funcPath = [targetF filesep Subjects{sj}];
                struct_dir = fullfile(funcPath, 'anat');
                data_dir = fullfile(funcPath, 'dec_SVR');
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
%                           Step 10: Smoothing
% =========================================================================

% smoothing --> prefix: s
if smoothing
    kernel_size = [3 3 3]; % smoothing kernel of 3mm

    for sj = 1:numel(Subjects)
        if ismember(sj, excludeSJ)
            continue;
        else
            display(['Step 10, Smoothing: ' Subjects{sj}])
            sj_dir = [targetF filesep Subjects{sj}];
            data_dir = [sj_dir filesep 'dec_SVR'];
            % A8_smoothing_run performs smoothing of functional image with
            % a pre-defined kernel using the spm_jobman with pre-defined parameters.
            A8_smoothing_run(data_dir, Subjects{sj}, ['^wres_zcorr.nii'], kernel_size);
        end
    end
end