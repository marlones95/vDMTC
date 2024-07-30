%% Bayesian Model Comparisons

clear
close all

single_subject = 1;

cd 'D:\vDMTC\Behavioural Analysis\Bayesian Model\VBA'
load('posteriors.mat'); % Posteriors and fitting stats of full model
close all

sjs = 36;

%%%% Single Subject model comparisons: %%%%%

% =========================================================================
%  Comparing the free energies of the "full model" to the “null-model” 
% =========================================================================

fitting_stats_full_model = fitting_stats;
    
if single_subject
    load('posteriors_null'); % load posteriors of the null model (which uses
    % Fechner adjusted physical differences between f1 and f2 to quantify
    % decisions)

    bayesianP = zeros(1,36); % vector containing Bayesian p-values
    BayesF    = zeros(1,36); % vector containing Bayes factors
    below_001 = 0;           % Count of all participants with p < 0.001
    below_01  = 0;           % Count of all participants with p < 0.01 & p >= 0.001
    below_05  = 0;           % Count of all participants with p < 0.05 & p >= 0.01
    not_sign  = [];          % Count of all participants with p >= 0.05

    no_evidence       = 0; % Number of participants with a Bayes Factor showing no evidence for full model
    weak_evidence     = 0; % Number of participants with a Bayes Factor showing weak evidence for full model    
    moderate_evidence = 0; % Number of participants with a Bayes Factor showing moderate evidence for full model
    strong_evidence   = 0; % Number of participants with a Bayes Factor showing strong evidence for full model

    for i = 1: sjs
        dF = fitting_stats_full_model(i).F - fitting_stats(i).F;
        bayesianP(i) = 1./(1+exp(dF));
        BayesF(i) = exp(dF);
        if bayesianP(i) < 0.001
            below_001 = below_001 +1;
        elseif bayesianP(i) < 0.01
            below_01 = below_01 + 1;
        elseif bayesianP(i) < 0.05
            below_05 = below_05 + 1;
        else
            not_sign = [not_sign;bayesianP(i)];
        end
        if BayesF(i) < 1
            no_evidence = no_evidence+1;
        elseif BayesF(i)>1 && BayesF(i)<3
            weak_evidence = weak_evidence+1;
        elseif BayesF(i)>3 && BayesF(i)<10
            moderate_evidence = moderate_evidence+1;
        else strong_evidence = strong_evidence+1;
        end
    end
end

% =========================================================================
%               random effects Bayesian model selection
% =========================================================================

% loading posteriors of the to-be-compared models
load('posteriors_bias_precision'); % biased simple model (estimating bias and likelihood variance)
fitting_stats_bias_precision = fitting_stats;
load('posteriors_precision');      % simple model (estimating only likelihood variance)
fitting_stats_precision = fitting_stats;
load('posteriors_precision_TOE');  % unbiased TOE model (likelihood variance and prior variance)
fitting_stats_precision_TOE = fitting_stats;

% Model selection
[posterior,out] = VBA_groupBMC([fitting_stats_full_model(:).F; fitting_stats_precision(:).F;...
    fitting_stats_bias_precision(:).F;fitting_stats_precision_TOE(:).F]);