% 
% 
% matlab_scripts.m
% These scripts run the misfit and uncertainty codes
% 
% 20170817 -- calvizuri <celso.alvizuri@unil.ch>

function[] = get_fmtu(evid, model, depth, wpol)
% update here
nsol_rand = 10000000;
nsol_grid = 22121190;
kval = 40; % max ~ 185
depstr = sprintf('%03d', depth);
ncell = 11; % default 11. NOTE (ncells, npts) = (11, 429);  (21, 1554)
input_psmecafix = sprintf('OUTPUT_DIR/%s_%s_%s_misfit_wf_psmeca', evid, model, depstr);

%c = parcluster
%set(c, 'NumWorkers', 16);
%-----------------------------------------------------------
% key routines
fmtu_post_process_misfit('OUTPUT_DIR', evid, model, depth, nsol_grid);
%psmeca_fix(input_psmecafix);   % 2019-08-02 has issues when eig=NAN or INF. appears an issue when too near lune boundaries (1/100 instead of 1/20)
fmtu_post_process_stp3('OUTPUT_DIR/', evid, model, depth, nsol_rand, kval, ncell, wpol);
fmtu_post_process_posterior(evid, 'OUTPUT_DIR/', model, depth, nsol_rand, kval, wpol, 'OUTPUT_DIR/');
%-----------------------------------------------------------
delete(gcp('nocreate'));

% rename output files
movefile(strcat('OUTPUT_DIR/', evid, '_histo2.dat')	          ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_histo2.dat'));                
movefile(strcat('OUTPUT_DIR/', evid, '_histo.dat')	          ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_histo.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_Mbeach_file.dat')	  ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_Mbeach_file.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_Mpost_file.dat')	      ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_Mpost_file.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_omega_err_post.dat')	  ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_omega_err_post.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_omega_misfit_outline.dat'), strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_omega_misfit_outline.dat'));  
movefile(strcat('OUTPUT_DIR/', evid, '_post_cdf.dat')	      ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_post_cdf.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_post_pdf.dat')	      ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_post_pdf.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_result.dat')	          ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_result.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_samples.dat')	      ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_samples.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_vP_avg.dat')	          ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_vP_avg.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_vP.dat')	              ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_vP.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_mesa_pdf.dat')	      ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_mesa_pdf.dat')); 
movefile(strcat('OUTPUT_DIR/', evid, '_mesa_cdf.dat')	      ,    strcat('OUTPUT_DIR/', evid, '_', model, '_', depstr, '_mesa_cdf.dat')); 
