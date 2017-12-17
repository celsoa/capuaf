function [v, w, kappa, theta, sigma, mag, misfit_wf, misfit_fmp, VR] = read_capbin_gd(inputfile, nrows)
% This script reads data output by CAP and loads it in MATLAB
% CAP output files are binary, the file name includes number of solutions (nrows)
%
%   capout_grid_bb_000100000.bin
%   capout_rand_mt_000100000.bin
%
% USAGE
%
% load data generated by CAP in GRID search
% [v, w, kappa, theta, sigma, mag, misfit_wf, misfit_fmp] = read_capbin_gd('capout_grid_gd.bin', nrows);
%
% load data generated by CAP in RANDOM search
% [v, w, kappa, theta, sigma, mag, misfit_wf, misfit_fmp] = read_capbin_gd('capout_rand_gd.bin', nrows);
%
% 20151216 celso alvizuri - cralvizuri@alaska.edu 
%----------------------------------------------------------

% read data
m = memmapfile(inputfile, 'Format', {'single', [11 nrows] 'floats'});
%m.data(1).floats(1:16);
a = m.data(1).floats()';
gamma      = a(:,1);
delta      = a(:,2);
kappa      = a(:,3);
theta      = a(:,4);
sigma      = a(:,5);
mag        = a(:,6);
misfit_fmp = a(:,7);
misfit_wf  = a(:,8);
        VR = a(:,9);
         v = a(:,10);
         w = a(:,11);
