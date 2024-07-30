function [cTOE,PerfGlobal] = Plot_TOE(u,y,flag)

% Script adapted from Jan Ostermann (born Herding)
%
% IN:
%   - u: design
%   - y: behavioral response
%   - flag: if 1 -> allow plot
% OUT:
%   - cTOE: slope and intercept of TOE fit
%   - PerfGlobal: global performance


% Get true frequency differences
F1F2 = u(1,:)-u(2,:);   % F1 - F2
IdF1F2 = find(F1F2>0);  % F1 > F2
IdF2F1 = find(F1F2<0);  % F1 < F2
IdRespF1F2 = find(y==1);% Response F1 > F2
IdRespF2F1 = find(y==0);% Response F1 < F2


% Get Accuracies for conditions
R = unique(u(1,:));
Nd = length(R);
TOE_Acc     = zeros(Nd,2);
for k = 1:Nd
    Istim = find(u(1,:) == R(k)); % Id F1
    I1 = intersect(Istim,IdF1F2);% Id F1>F2
    TOE_Acc(k,1)    = length(intersect(IdRespF1F2,I1))/length(I1);% Resp(F1>F2)
    I2 = intersect(Istim,IdF2F1);% Id F2>F1
    TOE_Acc(k,2)    = length(intersect(IdRespF2F1,I2))/length(I2);% Resp(F2>F1)
end

% Extract global Task Performance
PerfGlobal = (length(intersect(IdF1F2,IdRespF1F2))+length(intersect(IdF2F1,IdRespF2F1))) ...
    /length(F1F2);

% Plot accuracies
if flag
    % Display
    figure;
    set(gcf,'color','white');
    subplot(1,2,1);
    plot(R,TOE_Acc(:,1),'b-x','Linewidth',3);
    hold on
    plot(R,TOE_Acc(:,2),'r-.o','Linewidth',3);
    title('Accuracy','FontSize',20,'FontWeight','bold');
    axis([R(1) R(end) 0 1]);
    ylabel('% Correct','FontSize',20,'FontWeight','bold');
    xlabel('F1','FontSize',20,'FontWeight','bold');
    legend({'F2 < F1' 'F2 > F1'});
    set(gca,'FontSize',16,'FontWeight','bold');
else
end


% Get probability for choosing f1 > f2
R = unique(u','rows');
Nd = length(R);
TOE_Acc = zeros(Nd,1);
for k = 1:Nd
    Istim = intersect(find(u(1,:) == R(k,1)) , find(u(2,:) == R(k,2)));
    TOE_Acc(k) = length(find(y(Istim)==0))/length(Istim);
end

% Calculate bias for choosing f1 > f2
poss_f1 = unique(u(1,:));
Nf1     = length(poss_f1);
TOE_mean = zeros(Nf1,1);
for k =1:Nf1
    TOE_mean(k) = mean(TOE_Acc(poss_f1(k)==R(:,1)));
end

% Plot bias
if flag
    % Display
    subplot(1,2,2);
    plot(R(diff(R')<0,1),TOE_Acc(diff(R')<0),'bx','Linewidth',3);
    hold on
    plot(R(diff(R')>0,1),TOE_Acc(diff(R')>0),'ro','Linewidth',3);
    plot(poss_f1,TOE_mean,'k-s','Linewidth',3);
    title('Bias','FontSize',20,'FontWeight','bold');
    axis([R(1,1) R(end,1) 0 1]);
    ylabel('% resp(F2>F1)','FontSize',20,'FontWeight','bold');
    xlabel('F1','FontSize',20,'FontWeight','bold');
    legend({'F2 < F1' 'F2 > F1' 'F2=F1 estimation'});
    set(gca,'FontSize',16,'FontWeight','bold');
else
end

% Extract slope and intersect TOE
X = poss_f1;
Y = TOE_mean;
cTOE = robustfit(X,Y);

% Plot bias slope
if flag
    plot(X,cTOE(1)+cTOE(2)*X,'k','Linewidth',2);
else
end