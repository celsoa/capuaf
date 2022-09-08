function Mx = psmeca_fix_pertM(M)
%PSMECA_FIX_PERTM decrease fault dip angles by a small perturbation (GMT bug)
%
% For some moment tensors there are plotting bugs that arise that can be
% fixed by perturbing the moment tensor. Here we apply a perturbation by
% decreasing the dip angle of ALL moment tensors
% by a random number between PMIN and PMAX degrees
%
% The problem seems to arise for vertical faults.
%
% called by psmeca_fix.m

PMIN = 0.1;
PMAX = 0.2;

n = length(M);

% STEP 1: convert to FMT grid search parameters
[gamma,delta,M0,kappa,theta,sigma] = CMT2TT(M);

% STEP 2: decrease the dip angle
% (note: horizontal faults should not be allowed for this modification)
pvec = rand(n,1)*(PMAX-PMIN) + PMIN;
Mx = TT2CMT(gamma,delta,M0,kappa,theta-pvec,sigma);    % N-m
