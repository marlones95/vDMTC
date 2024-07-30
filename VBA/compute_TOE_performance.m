%% Compute d' and empirical Time Order Effect

function [global_dprime, bias_fit] = compute_TOE_performance(stim_pair, response)

% Script adapted from Jan Ostermann (born Herding)

% compute d' and the bias (criterion) of each subject for each value of f1
% define f1 > f2 trials as targets:
%
%                            stimulus pair
%                     f1>f2        |     f2>f1
%               ----------------------------------------
%         f1>f2 |     Hit (H)      |  false alarm (FA)
% choice  ------|---------------------------------------
%         f2>f1 |     Miss (M)     |  correct rejection emp_TOE(CR)
%               ----------------------------------------
%                      H + M       |      FA + CR
%
%   hit rate:          h = H/(H+M)
%   false alarm rate: fa = FA/(FA+CR)
%   d' = z(hit rate) - z(false alarm rate)
%
% -------------------------------------------------------------------------
%
% INPUT
%   stim_pair:      2-by-nTrials vector of stimulus frequencies
%   chose_f1:       1-by-nTrials binary vector with ones if response was f1>f2 and zeros else
%
% OUTPUT
%   global_dprime:  overall d'
%   bias_fit:       linear fit to the criterion value for each f1 -> measure of TOE!


% get the frequencies of f1 and the stim differences f2 - f1
f1         = stim_pair(1,:);
poss_f1    = unique(f1);
stim_diffs = diff(stim_pair);

% get trial indices with f1 > f2 & f2 > f1, respectively
f1_larger_f2 = find(stim_diffs < 0);
f2_larger_f1 = find(stim_diffs > 0);

% get trial indices of subject choices
chose_f1 = find(response==1);

% compute the performance for each f1 value
criterion = zeros(1,length(poss_f1));


for f1_idx=1:length(poss_f1)

    % get all trials with current f1 value
    curr_f1_trials = find(f1 == poss_f1(f1_idx));

    % get all the trials of current f1 value in which f1>f2 was true ...
    f1_larger_curr_f1 = intersect(curr_f1_trials, f1_larger_f2);
    % ... and in which f2>f1 was true
    f2_larger_curr_f1 = intersect(curr_f1_trials, f2_larger_f1);

    % hit_rate
    h = numel(intersect(f1_larger_curr_f1, chose_f1))/...  % number of hits for given f1 value
        numel(f1_larger_curr_f1);                          % number of trials with f1>f2 for given f1

    % adjust hit rate if it was perfect
    if h == 1
        h = 1 - 1/(2*numel(f1_larger_curr_f1));
    end

    % FA rate
    fa = numel(intersect(f2_larger_curr_f1, chose_f1))/...  % number of false alarms (FA) for given f1
        numel(f2_larger_curr_f1);                           % number of trials with f2>f1 for given f1

    % adjust false alarm rate if it was perfect
    if fa == 0
        fa = 1/(2*numel(f2_larger_curr_f1));
    end

    % criterion c = -0.5*(z(h)+z(fa))
    criterion(f1_idx) = -mean([norminv(h) norminv(fa)]);

end

% hit rate
global_h = numel(intersect(f1_larger_f2, chose_f1))/...     % number of hits
    numel(f1_larger_f2);                                    % number of trials with f1>f2

% FA rate
global_fa = numel(intersect(f2_larger_f1, chose_f1))/...    % number of false alarms (FA)
    numel(f2_larger_f1);                                    % number of trials with f2>f1

% d' measure
global_dprime = (norminv(global_h) - norminv(global_fa))/sqrt(2);

% Empirical TOE
bias_fit = robustfit(poss_f1,criterion(1,:));
end