%% Observation Function
function [gx] = g_vDMTC_simple_bias(x,Phi,~,~)

% Script adapted from Jan Ostermann (born Herding)
%
% Computes the probability of the subject choosing f1 > f2 and
% derives the action of the subject following a softmax rule.
% IN:
%   - x: the state posterior sufficient statistics.
%   - Phi: Phi-parameters
% OUT:
%   - gx: the predicted subject's output.

% up-dated conditional sufficient statistics on states
m1p = x(1);
v1p = x(2);
m2p = x(3);
v2p = x(4);
m1  = x(5);
v1  = x(6);
m2  = x(7);
v2  = x(8);

bias = Phi(1);

% =========================================================================
%                       Compute choice probability
% =========================================================================

Dm = m1p - m2p;     % difference between posterior means
Sv = v1p + v2p;     % sum of posterior variances


gx = 0.5*(1 + erf(Dm/(sqrt(2)*sqrt(Sv)))) + bias; % Prob F1>F2

% Y = ERF(X) is the error function defined as:
% erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.

end