%% Evolution Function
function [fx] = f_vDMTC_simple(x,P,u,in)

% Script adapted from Jan Ostermann (born Herding)
%
% This is the update rule for the posterior sufficient statistics of a delayed match-to-comparison (DMTC) task.
%
% IN:
%   - x: states
%       m1 = x(1);
%       v1 = x(2);
%       m2 = x(3);
%       v2 = x(4);
%       m1 = x(5);
%       v1 = x(6);
%       m2 = x(7);
%       v2 = x(8); 
%   - P: the perceptual model parameter vector (in log space), ie. P = [sigma,mu_global,v_global]
%   - u: the two trial inputs (frequency of stim 1 and stim 2 in log space)
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u),
%         as well as book keeping of the past believes.


% =========================================================================
%                           Initialize values
% =========================================================================

% stimulus vals
f1 = u(1);          % f1 (freq on log scale)
f2 = u(2);          % f2 (freq on log scale)

% model parameters
sigma     = P(1);      % precision of encoding (variance of stimulus likelihood for f1 and f2)
mu_global = P(2);      % mean of global prior (on log scale)
v_global  = P(3);      % variance of global prior

% make sure that prior variance is interpretably (v_global>0)
if v_global <= 0
    v_global = nan;
end

% model states
m1  = mu_global;
v1  = v_global;
v_lhood1 = sigma;

% =========================================================================
%                           Update hidden states
% =========================================================================

%%%% Updating the hidden states pertaining to the first stimulation %%%


v1p = 1/(1/v1 + 1/v_lhood1);         % variance of posterior for f1
m1p = v1p * (m1/v1 + f1/v_lhood1);   % mean of posterior for f1 (log scale)


%%% Updating the hidden states pertaining to the second stimulation %%%

v_lhood2 = sigma;   % variance of stimulus likelihood (same for f1 and f2)

switch in.PriorF2
    case 'same'
        m2  = m1;
        v2  = v1;
        v2p = 1/(1/v2 + 1/v_lhood2);        % variance of posterior for f2
        m2p = v2p * (m2/v2 + f2/v_lhood2);  % mean of posterior for f2 (freq on log scale)
    case 'PostF1'
        m2  = m1p;
        v2  = v1p;
        v2p = 1/(1/v2 + 1/v_lhood2);        % variance of posterior for f2
        m2p = v2p * (m2/v2 + f2/v_lhood2);  % mean of posterior for f2 (freq on log scale)
    case 'none'
        m2  = m1;                           % set to arbitrary value
        v2  = v1;                           % set to infinty for display purposes only... prior not used
        m2p  = f2;                          % mean of likelihood of f2 serves as mean of posterior for f2
        v2p  = v_lhood2;                    % variance likelihood of f2 serves as variance of posterior for f2
end

%%% Save Conditional posteriors on hidden states in fx structure %%%
fx = zeros(8,1);
fx(1) = m1p;
fx(2) = v1p;
fx(3) = m2p;
fx(4) = v2p;
fx(5) = m1;
fx(6) = v1;
fx(7) = m2;
fx(8) = v2;

end