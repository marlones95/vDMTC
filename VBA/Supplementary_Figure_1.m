%% Supplementary Figure 1

% In this script a plot with 36 subplots is made. Each plot compares the
% prediction of the model (squares) to the true data (lines) for f1 < f2
% (blue) and f1 > f2 (red). The x-axis represents f1 magnitude, the y-axis
% performance

% Required script: tight_subplot.m 
% https://de.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w

clear
close all

% Subjects
SJs = {'01', '02', '03', '04', '05','06', '07', '08', '09', '10','11', '12',...
    '13',  '14','15', '16', '17', '18','19', '20','21', '22', '23','24', '26', '27',...
    '28','29', '30','31', '32', '33', '34', '35', '36', '37'};

% log file directory
data_dir = 'D:\vDMTC\Logs';
sjs = size(SJs,2);

% Load posteriors
load('posteriors.mat');

stim_set = [16 16 20 20 24 24 28 28;...
    12 20 16 24 20 28 24 32];
log_stim_set = log(stim_set);

f1 = log_stim_set(1,:);          % f1 (freq on log scale)
f2 = log_stim_set(2,:);          % f2 (freq on log scale)

% mean of global prior (on log scale)
mu_global = mean(log_stim_set(:));  % prior mean


figure
% Create subplots with reduced spacing
gap = [0.02 0.02];      % Gap between subplots [horizontal, vertical]
marg_h = [0.1 0.05];    % Margins [bottom, top]
marg_w = [0.05 0.01];   % Margins [left, right]
ha = tight_subplot(6, 6, gap, marg_h, marg_w);


% Dependent Variable: Performance per condition
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


    perf_F1_F2 = [perf_16_12, perf_20_16, perf_24_20, perf_28_24];
    perf_F2_F1 = [perf_16_20, perf_20_24, perf_24_28, perf_28_32];

    % =====================================================================
    %           Estimate predicted performance per condition
    % =====================================================================
    
    m1       = mu_global; % prior mean
    v1       = subj_posterior(i).muTheta(4); % prior variance
    v_lhood1 = subj_posterior(i).muTheta(1); % likelihood variance of f1
    v_lhood2 = subj_posterior(i).muTheta(1); % likelihood variance of f2
    v2p      = v_lhood2;                     % posterior variance of f2 (equals likelihood variance)
    bias     = subj_posterior(i).muPhi;      % response bias

    gx = zeros(1,8); % Prob F1>F2 for all conditions                       

    for j = 1:length(stim_set)
        v1p = 1/(1/v1 + 1/v_lhood1);            % variance of posterior for f1
        m1p = v1p * (m1/v1 + f1(j)/v_lhood1);   % mean of posterior for f1 (freq on log scale)
        m2p = f2(j);                            % mean posterior of f2 (equals likelihood mean)

        %%% Estimate behavioural responses %%%

        Dm       = m1p - m2p;                             % difference of posterior means
        Sv       = v1p + v2p;                             % sum of posterior variance
        chose_f1 = 0.5*(1 + erf(Dm/(sqrt(2)*sqrt(Sv))));  % Prob F1>F2
        gx(j)    = chose_f1 + bias;                       % Prob F1>F2 including bias

        subj_diff = exp(m2p)-exp(m1p); % subjective difference in Hz

    end

    % Caclulate estimated performances
    perf_F1_F2_model = [gx(1),gx(3),gx(5),gx(7)];
    perf_F2_F1_model = [1-gx(2),1-gx(4),1-gx(6),1-gx(8)];

    xCnt = [0.05, 0.35, 0.65, 0.95];

    % =====================================================================
    %               Plot predicted and true performances
    % =====================================================================
    
    % Plot true performances
    axes(ha(i));
    hold on
    p(1) = plot(xCnt,perf_F1_F2,'rs','MarkerSize',3);
    p(2) = plot(xCnt,perf_F2_F1,'bs','MarkerSize',3);
    set(p(1), 'markerfacecolor', get(p(1), 'color'));
    set(p(2), 'markerfacecolor', get(p(2), 'color'));

    % Plot predicted performances
    p(3) = plot(xCnt,perf_F1_F2_model,'r-','LineWidth',2);
    p(4) = plot(xCnt,perf_F2_F1_model,'b-','LineWidth',2);
    xlim([0,1]);
    ylim([0.25,1.1]);
    xticks(sort(xCnt(:)))
    yticks([0.25 0.5,0.75,1]);
    set(gca,'linewidth',1)
    hold off
end



