%% Figure 2D

clear all
close all


SJs = {'01', '02', '03', '04', '05','06', '07', '08', '09', '10','11', '12',...
    '13',  '14','15', '16', '17', '18','19', '20','21', '22', '23','24', '26', '27',...
    '28','29', '30','31', '32', '33', '34', '35', '36', '37'};
data_dir = 'D:\vDMTC\Logs'

sjs = size(SJs,2);


% Load posteriors from the VBA
cd("D:\vDMTC\Behavioural Analysis\Bayesian Model\VBA")
load('posteriors.mat');

stim_set = [16 16 20 20 24 24 28 28;...
    12 20 16 24 20 28 24 32];
log_stim_set = log(stim_set);

f1 = log_stim_set(1,:); % f1 (freq on log scale)
f2 = log_stim_set(2,:); % f2 (freq on log scale)

% Setting up matrix for true performances
F1_higher = zeros(sjs, 4);  % true performance of conditions with f1>f2
F1_lower  = zeros(sjs, 4);  % true performance of conditions with f1<f2

% Bayesian Parameters
sigma       = zeros(1,sjs);          % likelihood variance
mu_global   = mean(log_stim_set(:)); % mean of global prior (on log scale)
v_global    = zeros(1,sjs);          % variance of global prior
bias        = zeros(1,sjs);          % response bias

% Setting up matrix for predicted performances
gx_all = zeros(2,4,36);

for i = 1:sjs
    sj = SJs{i};
    logfile = spm_select('FPList', data_dir, ['^fMRI_', sj, '.*\.mat$']);
    load(logfile);
    % =====================================================================
    %               Calculate performance for each condition
    % =====================================================================

    perf_16_12 = sum(mylog.perform(mylog.flutter(:,:,1) == 16 & mylog.flutter(:,:,2) == 12))/48;
    perf_16_20 = sum(mylog.perform(mylog.flutter(:,:,1) == 16 & mylog.flutter(:,:,2) == 20))/48;
    perf_20_16 = sum(mylog.perform(mylog.flutter(:,:,1) == 20 & mylog.flutter(:,:,2) == 16))/48;
    perf_20_24 = sum(mylog.perform(mylog.flutter(:,:,1) == 20 & mylog.flutter(:,:,2) == 24))/48;
    perf_24_20 = sum(mylog.perform(mylog.flutter(:,:,1) == 24 & mylog.flutter(:,:,2) == 20))/48;
    perf_24_28 = sum(mylog.perform(mylog.flutter(:,:,1) == 24 & mylog.flutter(:,:,2) == 28))/48;
    perf_28_24 = sum(mylog.perform(mylog.flutter(:,:,1) == 28 & mylog.flutter(:,:,2) == 24))/48;
    perf_28_32 = sum(mylog.perform(mylog.flutter(:,:,1) == 28 & mylog.flutter(:,:,2) == 32))/48;


    F1_higher(i,:) = [perf_16_12, perf_20_16, perf_24_20, perf_28_24];
    F1_lower(i,:)  = [perf_16_20, perf_20_24, perf_24_28, perf_28_32];

    % =====================================================================
    %           Estimate predicted performance per condition
    % =====================================================================
    
    % model parameters
    sigma(i)     = subj_posterior(i).muTheta(1);
    v_global(i)  = subj_posterior(i).muTheta(4);
    bias(i)      = subj_posterior(i).muPhi;

    m1       = mu_global; % prior mean
    v1       = subj_posterior(i).muTheta(4); % prior variance
    v_lhood1 = subj_posterior(i).muTheta(1); % likelihood variance of f1
    v_lhood2 = subj_posterior(i).muTheta(1); % likelihood variance of f2
    v2p      = v_lhood2;                     % posterior variance of f2 (equals likelihood variance)
    b        = subj_posterior(i).muPhi;      % response bias

    gx = zeros(1,8); % Prob F1>F2 for all conditions                       

    for j = 1:length(stim_set)
        v1p = 1/(1/v1 + 1/v_lhood1);            % variance of posterior for f1
        m1p = v1p * (m1/v1 + f1(j)/v_lhood1);   % mean of posterior for f1 (freq on log scale)
        m2p = f2(j);                            % mean posterior of f2 (equals likelihood mean)

        %%% Estimate behavioural responses %%%

        Dm       = m1p - m2p;                             % difference of posterior means
        Sv       = v1p + v2p;                             % sum of posterior variance
        chose_f1 = 0.5*(1 + erf(Dm/(sqrt(2)*sqrt(Sv))));  % Prob F1>F2
        gx(j)    = chose_f1 + b;                       % Prob F1>F2 including bias

        subj_diff = exp(m2p)-exp(m1p); % subjective difference in Hz

    end

    % Caclulate Prob F1>F2
    gx_all(:,:,i) = [gx(1),gx(3),gx(5),gx(7);...
        1-gx(2),1-gx(4),1-gx(6),1-gx(8)];
end

means_F1_higher = mean(F1_higher);
means_F1_lower  = mean(F1_lower);

means  = mean(gx_all,3);
means_f1_higher_m = [means(1,:)];
means_f1_lower_m  = [means(2,:)];




sterrs = std(gx_all,[],3)/sqrt(size(gx_all,3));
sterrs = [sterrs(1,:),sterrs(2,:)];

xCnt = [0.05, 0.35, 0.65, 0.95];

% =========================================================================
%                   Plot predicted and true performances
% =========================================================================

figure
hold on

% Generate dummy handles for legend
h(1) = plot(xCnt,means_F1_higher,'ks');
set(h(1), 'markerfacecolor', get(h(1), 'color'));
h(2) = plot(xCnt,means_f1_lower_m,'k-','LineWidth',3);
h(3) = plot(xCnt,means_F1_higher,'rs');
set(h(3), 'markerfacecolor', get(h(3), 'color'));
h(4) = plot(xCnt,means_F1_lower,'bs');
set(h(4), 'markerfacecolor', get(h(4), 'color'));

% Plot true performances
p(1) = plot(xCnt,means_F1_higher,'rs');
p(2) = plot(xCnt,means_F1_lower,'bs');
set(p(1), 'markerfacecolor', get(p(1), 'color'));
set(p(2), 'markerfacecolor', get(p(2), 'color'));

% Plot predicted performances
p(3) = plot(xCnt,means_f1_higher_m,'r-','LineWidth',4);
p(4) = plot(xCnt,means_f1_lower_m,'b-','LineWidth',4);
xlim([0,1]);
ylim([0.5,1]);
xticks(sort(xCnt(:)))
xticklabels({'16', '20', '24', '28','16', '20', '24', '28'});
yticks([0.5,0.75,1]);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',25)
set(gca,'linewidth',2)
hold off
legend(h,'data','model','f1 > f2','f1 < f2','Location','southwest','FontSize',16)