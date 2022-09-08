function M = psmeca_fix_increaseM(M)
%PSMECA_FIX_INCREASEM increase moment tensor magnitude to overcome GMT bug
%
% called by psmeca_fix.m

% increase ALL magnitudes to avoid the small-magnitude bug
IEXP = 10;
warning('make sure that psmeca_legend.m has the same IEXP = %i',IEXP);
M = M * 10^IEXP;
