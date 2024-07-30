%% Figure 2B and 2C
% Create 2 figures
% Figure 1: Correlates d' with likelihood precision
% Figure 2: Correlates empirical Time Order Effect with prior precision
clear
close all

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17',...
    '18','19','20','21','22','23','24','26','27','28','29','30','31','32','33','34','35','36','37'};

N = length(subjects);

log_dir = 'D:\vDMTC\Logs';

filenames = cell(1,length(subjects));
for i=1:length(subjects)
    filenames{i} = spm_select('FPList', log_dir, append('^fMRI_', subjects(i), '.*\.mat$'));
end


%%
% stimulus params
%------------

% load posteriors
load('posteriors.mat');
close all

% =========================================================================
%       Loop over subjects to obtain d' and empirical TOE
% =========================================================================

global_dprime = zeros(1, N);    % d'
bias_fit      = zeros(2, N);    % empirical TOE

for n = 1:N
    
    % empirical data
    load(fullfile(filenames{n}))

    stim_pairs      = reshape(permute(mylog.flutter, [3 2 1]), 2, []);

    f1_f2_correct   = mylog.perform == 1 & (mylog.flutter(:,:,1) > mylog.flutter(:,:,2));
    f1_f2_incorrect = mylog.perform == 0 & (mylog.flutter(:,:,1) < mylog.flutter(:,:,2));
    mylog.resp      = f1_f2_incorrect + f1_f2_correct;
    mylog.resp(isnan(mylog.choice_content)) = -1;
    behav_resp      = reshape(mylog.resp', 1, []); % 1=f1>f2 & 0=f2>f1 & -1 = no answer

    valid_idx = all([stim_pairs(1,:) ~= 0 ; behav_resp ~= -1]);
   
    emp_u = stim_pairs(:,valid_idx);    
    emp_y = behav_resp(valid_idx);   
    
    % Get measure of d' and emprical TOE
    [global_dprime(n),bias_fit(:,n)] = compute_TOE_performance(emp_u,emp_y);
end

% =========================================================================
%                           Make figures
% =========================================================================

%%% Figure 2b: Correlate d' with estimated likelihod precision %%%
f2_col =[0 0.8 1];
thetas     = [subj_posterior.muTheta];
precisions = 1./thetas(1,:);
prior_precision = 1./thetas(3,:);

% Make figure
figure
scatter(global_dprime,precisions,30,f2_col,'filled')
[R,P] = corrcoef(global_dprime,precisions);
h = lsline;
h.Color = f2_col;
h.LineWidth = 4;
set(gca,'fontname','arial')
xticks([1,1.4,1.8,2.2,2.6]);
yticks([0,100,200]);
xlim([1,2.6]);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',25)
set(gca,'linewidth',2)

%%% Figure 2c: Correlate empirical TOE with prior precision (true TOEvs. estimated TOE) %%%
emp_TOE = bias_fit(2,:);

% Make figure
figure
scatter(emp_TOE,prior_precision,30,'k','filled')
[R_prior,P_prior] = corrcoef(emp_TOE,prior_precision); 
h = lsline;
h.Color = 'k';
h.LineWidth = 4;
set(gca,'fontname','arial')
xticks([-0.04,0,0.04,0.08,0.12]);
yticks([0,40,80]);
xlim([-0.04,0.12]);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',25)
set(gca,'linewidth',2)