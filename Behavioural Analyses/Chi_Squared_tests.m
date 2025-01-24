%% Chi Squared Tests:
% This script performs Pearson chi-square tests comparing the proportions of the two
% conditions "higher" vs "lower" between left and right motor responses
clear

% Required script: prop_test.m
% https://de.mathworks.com/matlabcentral/fileexchange/45966-compare-two-proportions-chi-square

% Initialize Subject cell array
SJs     = {'sub-002','sub-003','sub-004','sub-005','sub-007','sub-008','sub-009','sub-010','sub-012','sub-016','sub-017',...
    'sub-018','sub-022','sub-023','sub-026','sub-027','sub-029','sub-030','sub-031','sub-032','sub-033','sub-034','sub-036','sub-037'};
% SJ25 excluded because of problems with projector
% Subjects excluded because performance in one condition < 0.5: SJ01, SJ06, SJ11, SJ13, SJ14, SJ19, SJ21, SJ24, SJ28
% Subjects excluded with head movement > 3mm: SJ15, SJ20, SJ35

% Directory containing log files
log_dir = 'D:\vDMTC\Logs';

% Initialize relevant vectors
h_motor = zeros(length(SJs),1); % hypothesis vector (0 = distributions not significantly different,
% 1 = distributions significantly different)
p_motor = zeros(length(SJs),1); % vector containing the p-value of the chi-squared test
chi2stat_motor = zeros(length(SJs),1); % vector containing the chi-squared statistic

for sj = 1:length(SJs)
    current_sj = SJs{sj}(end-1:end);
    logfile = spm_select('FPList', log_dir, ['^fMRI_', current_sj, '.*\.mat$']);
    load(logfile);

    decision_high_left = (mylog.choice_direction == 1 & ((mylog.focus == 1 & mylog.flutter(:,:,1) >...
        mylog.flutter(:,:,2)) | (mylog.focus == 2 & mylog.flutter(:,:,1) < mylog.flutter(:,:,2))));
    decision_high_right = (mylog.choice_direction == 2 & ((mylog.focus == 1 & mylog.flutter(:,:,1) >...
        mylog.flutter(:,:,2)) | (mylog.focus == 2 & mylog.flutter(:,:,1) < mylog.flutter(:,:,2))));
    decision_low_left = (mylog.choice_direction == 1 & ((mylog.focus == 2 & mylog.flutter(:,:,1) >...
        mylog.flutter(:,:,2)) | (mylog.focus == 1 & mylog.flutter(:,:,1) < mylog.flutter(:,:,2))));
    decision_low_right = (mylog.choice_direction == 2 & ((mylog.focus == 2 & mylog.flutter(:,:,1) >...
        mylog.flutter(:,:,2)) | (mylog.focus == 1 & mylog.flutter(:,:,1) < mylog.flutter(:,:,2)))); 

    X_left = [sum(sum(decision_high_left)),sum(sum(decision_low_left))];
    X_right = [sum(sum(decision_high_right)),sum(sum(decision_low_right))]; 
    N = [X_left(1)+X_right(1), X_left(2)+X_right(2)];

    % Chi Squared Tests
    [h_motor(sj),p_motor(sj), chi2stat_motor(sj)] = prop_test(X_left,N,'true');

end

% Creating a table containing all outputs of the chi-squared tests per
% participant
T_motor = table(SJs',h_motor,p_motor,chi2stat_motor);
disp(T_motor)