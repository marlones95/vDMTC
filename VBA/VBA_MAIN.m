%% VBA_model

% Script adapted from Jan Ostermann (born Herding)

%--------------------------------------------------------------------------
% This script models perceptual decisions in a  delayed match-to-comparison
% (DMTC) task.
% Task : 2 sequentially presented flicker stimuli (here visual flicker)
% have to be compared according their frequency
% The model accounts for learning from past observations and the time order
% effect which arises from maintaining a mnemonic representation of a
% stimulus throughout a delay-period. The Model updates the belief based on
% perceived stimuli and predict a response at the same trial.
%--------------------------------------------------------------------------

clear
close all

% =========================================================================
%       Initialize Subjects, Directories, and Parameters for Analysis
% =========================================================================

% list of subjects
subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13',...
    '14','15','16','17','18','19','20','21','22','23','24','26','27','28','29',...
    '30','31','32','33','34','35','36','37'};

% Log file directory
log_dir = 'D:\vDMTC\Logs';

% Extract filenames
filenames = cell(1,length(subjects));
for i=1:length(subjects)
    filenames{i} = spm_select('FPList', log_dir, append('^fMRI_', subjects(i), '.*\.mat$'));
end


%%% stimulus parameters %%%
% Initialize stimulus set
stim_set = [16 16 20 20 24 24 28 28;... % Frequency of first stimulus (f1)
    12 20 16 24 20 28 24 32];   % Frequency of second stimulus (f2)

% log transform the frequencies to account for Fechner's law
log_stim_set = log(stim_set);
% Get the mean of the log-transformed frequencies
f_mean   = mean(log_stim_set(:));


save_posteriors = 1; % Save posteriors in a .mat file
plot_results    = 0; 


% =========================================================================
%                           parameters & priors
% =========================================================================

%% Initialize priors of estimated parameter
% likelihood variance
sigma    = 0.05;    % Initial variance parameter of likelihood (default: 0.05)
mu_prior = f_mean;  % Initialize the mean of the global prior as mean over log stimulus set

% Initial variance of global prior: Full Width at Half Maximum of the range of all frequencies
fwhm         = max(log_stim_set(2,:))-min(log_stim_set(2,:));
v_prior      = (fwhm/(2*sqrt(2*log(2))))^2;

% Initial bias term
bias = 0;

% Evolution Parameters (1 = estimated, 0 = fixed)
param2est     = [1 0 1];    % sigma, mu_prior, v_prior
% Observation Parameters (1 = estimated, 0 = fixed)
phi_param2est = [1];        % bias

est_var       = 1;          % variance of to-be-estimated posterior distributions
param2est     = est_var.*param2est;
phi_param2est = est_var.*phi_param2est;

%% Evolution and Observation functions
f_fname = @f_vDMTC_simple;       % evolution function perception process
g_fname = @g_vDMTC_simple_bias;  % whole decision process


criterion_compare = zeros(5,2,27);

% =========================================================================
%               Loop over all subjects & model inversion
% =========================================================================

% Loop over all subjects
for n = 1:length(filenames)

    load(filenames{n});

    % Reshape frequency pairs, adding together all runs to a 2*384matrix
    stim_pairs      = reshape(permute(mylog.flutter, [3 2 1]), 2, []);
    % Extract responses (1 = f1 > f2, 0 = f1 < f2, nan = no response)
    behav_resp_all  = NaN(size(mylog.perform));
    behav_resp_all(mylog.perform == 1 & mylog.flutter(:,:,1)>mylog.flutter(:,:,2)) = 1;
    behav_resp_all(mylog.perform == 0 & mylog.flutter(:,:,1)<mylog.flutter(:,:,2)) = 1;
    behav_resp_all(mylog.perform == 1 & mylog.flutter(:,:,1)<mylog.flutter(:,:,2)) = 0;
    behav_resp_all(mylog.perform == 0 & mylog.flutter(:,:,1)>mylog.flutter(:,:,2)) = 0;
    behav_resp_all(isnan(mylog.choice_content)) = nan;
    behav_resp_all  = reshape(behav_resp_all',1, []);
    % Erase all non-responses
    behav_resp      = behav_resp_all(~isnan(behav_resp_all));
    stim_pairs      = stim_pairs(:,~isnan(behav_resp_all));

    u = stim_pairs; % design (the data)
    y = behav_resp; % the input (behaviour)

    %%% prior mean set %%%
    m1 = f_mean;    % prior mean for f1
    m2 = m1;        % prior mean for f2
    v1 = v_prior;   % prior variance f1
    v2 = v1;        % prior variance f2

    % Option for f function
    inF.PriorF2   = 'none'; % 'none' (use likelihood of f2), 'PostF1' (use posterior of f1) or
    % 'Same' (use same as for f1)

    % =====================================================================
    %           Initialization of parameters (Theta, Phi, x0)
    % =====================================================================
    theta = [sigma; mu_prior; v_prior]; % f parameters: sigma and v_prior are estimated
    phi   = bias; % g parameters: bias is estimated
    % Print intial values
    fprintf('sigma:\t%1.3f\nmu_p:\t%1.3f\nv_p:\t%1.3f\nbias:\t%1.3f\n', theta, phi)

    x0 = repmat([m1;v1;m2;v2],2,1); % hidden states priors

    % choose initial conditions (priors)
    dim = struct('n',length(x0),...     % number of hidden states in the probability learning model
        'n_theta',length(theta),...     % evolution parameters
        'n_phi',numel(phi));            % observation parameters

    % Assign values to options input structure
    options = [];
    options.sources.type = 1;
    priors.muPhi         = phi;
    priors.muTheta       = theta;
    priors.muX0          = x0;
    priors.SigmaX0       = 0*eye(dim.n);
    priors.SigmaPhi      = 0*eye(dim.n_phi);
    priors.SigmaTheta    = 0*eye(dim.n_theta);
    priors.a_alpha       = Inf;
    priors.b_alpha       = 0;

    options.priors   = priors;
    options.inF      = inF;

    options.priors.a_alpha  = Inf;
    options.priors.b_alpha  = 0;
    options.priors.muX0     = x0;
    options.priors.SigmaX0  = diag([0 5 0 5 0 5 0 5]);

    options.priors.muPhi    = phi;
    options.priors.SigmaPhi = diag(phi_param2est);

    options.priors.muTheta    = theta;
    options.priors.SigmaTheta = diag(param2est);

    options.isYout = zeros(size(y,1),size(y,2));

    options.backwardLag = 1; % size of short sighted backward passes (default: 1)
    %----------------------------------------------------------------------
    %                           MODEL INVERSION
    %----------------------------------------------------------------------
    [posterior,out] = VBA_NLStateSpaceModel(y,log(u),f_fname,g_fname,dim,options);

    %%% Plot results if plot_results = 1 %%%
    [cTOE(n,:),PerfGlobal] = Plot_TOE(u,y,plot_results); % Plot Accuracies and bias
    
    % display posterior means of theta parameters  
    if plot_results
        
        xtick = 1:out.dim.n_theta;
        if ~out.options.OnLine
            V = VBA_getVar(posterior.SigmaTheta);
            muTheta = posterior.muTheta;
        else
            V = VBA_getVar(posterior.SigmaTheta{end});
            muTheta = posterior.muTheta(:,end);
        end
        plotUncertainTimeSeries(muTheta,V,[]);
        set(gca, 'Xtick', xtick)
        title('theta')
        plot(theta,'go')
    end

    %%% Save model inversion results %%%
    subj_posterior(n) = posterior;
    fitting_stats(n)  = out;
    subj_ID(n)        = str2num(subjects{n});

end

% Save posteriors
if save_posteriors
    close all
    save("posteriors.mat","subj_posterior","fitting_stats")
end