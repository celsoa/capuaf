#!/usr/bin/env perl
# if on pacman uncomment the following
#use v5.10;
#
# A user-friendly PERL interface to the CAP source inversion code cap
#
# written by Lupei Zhu, 3/6/1998, Caltech
# 
# revision history
#	6/18/2001	add usage and documentation.
#	11/05/2012	add isotropic and CLVD search (-J).
#

use File::Copy;

# these are the only things one need to change based on the site installation
$home = $ENV{HOME};                        # my home directory
$caphome = $ENV{CAPHOME};                  # CAP home directory
$caprun = $ENV{CAPRUN};                    # run directory

require "$caphome/cap_plt.pl";             # include plot script
require "$caphome/sub_read_parameter_file.pl";  #  read parameter file

# rename output dir to include time. Options: 
# 1. propagate into CAP.c. 
# 2: simply rename capc's output_dir into this dir.
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
$runtime1 = sprintf("%4d-%02d-%02d", $year+1900, $mon, $mday);
$runtime2 = sprintf("%4d-%02d-%02dT%02d.%02d.%02d", $year+1900, $mon, $mday, $hour, $min, $sec);

#================defaults======================================
$cmd = "cap";

$inp_cmd = "inp_cmd";

# green's function location
# within each of these directories are subdirectories of models (cus, scak, utuhalf, etc)
#$green = "$home/data/models/Glib";               # original
$green = "/store/wf/FK_synthetics";               # UAF linux network
#$green = "/import/c/d/ERTHQUAK/FK_synthetics ";  # UAF cluster
#$green = "$caprun/models";                       # user testing
#$green = "/home/alvizuri/PROJECTS/CAP/models";   # western US model (temporary location: move to /store)
#$green = "/home/vipul/REPOSITORIES/capuaf/divider_things";

$repeat = 0;           # repeat inversion and discard bad trace (OBSOLETE)
$fm_thr = 0.01;        # first motion threshold
$dirct='';             # rupture directivity
$rise = 0.5;
$disp=0;               # integrate (velocity to displacment)
$pol_wt=999;           # relative weight to polarity misfit
$weight="weight.dat";  # deafult name for weight file
$fmt_flag="false";     # use 'fmt' flag for full moment tensor search outputs

# plotting
$plot = 0;             # to generate plots (waveform misfit, beachballs)
$ampbody_input = 1;        # amplitude scaling for P wave
$ampsurf_input = 2;        # amplitude scaling for Surf wave
$spib = 40;            # sec per inch, body waves
$spis = 45;            # sec per inch, surface waves
$keep = 0;             # keep synthetics

# filters and window lengths
($f1_pnl, $f2_pnl, $f1_sw, $f2_sw, $tmax_body, $tmax_surf) = (0.02,0.2,0.02,0.1,35,70);

# max. shifts
$max_shft1=1;		# max. shift for Pnl
$max_shft2=5;		# max. shift for surface wave
$Sstatic_shift = 0;     # default static shift for surface waves
$tie = 0.5;		# tie between SV and SH

# weights between different portions
$weight_of_pnl=2.0;	# weight for pnl window
$weight_of_rayleigh=1.0;# weight for rayeligh window
$weight_of_love=1.0;    # weight for love window
$power_of_body=1;	# distance scaling power for pnl waves
$power_of_surf=0.5;     # distance scaling power for surf waves

# apparent velocities
#($vp, $love, $rayleigh) = (7.8, 3.5, 3.1);
($vp, $love, $rayleigh) = (-1, -1, -1);

# search types (Random or grid)
$grid_type = -1.0;

# minimization (norm)
$norm = 1;

# for sorting the output file by distance or azimuth
$isort = 0;

# Use parameter file instead of command line input
$parameter_file = '';
$capver = ""; # default uses "cap"
#----------------------------------------------------------- 
# DEFAULT VALUES FOR SEARCH PARAMETERS
use Math::Trig ':pi';
$deg2rad = pi / 180.0;

# RANGE LIMITS PER PARAMETER
# magnitude 
($deg, $dm, $dlune) = (10, 0.1, 0.);
$nmw = 1;

# orientation
$k1 = 0;    $k2 = 360;
$h1 = 0;    $h2 = 1;
$s1 = -90;  $s2 = 90;

# lune --> (v,w) 
# v = [-1/3, 1/3] w = [-3pi/8, 3pi/8]
$w2 = (3.0 * pi / 8.0); $w1 = -$w2;
$v2 = (1.0 / 3.0);      $v1 = -$v2;

# NUMBER OF POINTS PER PARAMETER

# search type = GRID -- specify number of points for each parameter
# case 1. full moment tensor. use nv/nw/nstr/ndip/nrak
# case 2. fixed lambda. use nstr/ndip/nrak
# default number of points for each parameter

# for nX, X>1 (at least one parameter value has to exist!)
$nv = 1; $dv = 10;    # gamma
$nw = 1; $dw = 10;    # delta
$nk = 1; $dk = 10;    # strike
$nh = 1; $dh = 10;    # dip
$ns = 1; $ds = 10;    # rake

# magnitude default is to run a single point
$nmw = 1;
$dmw = 0;

# search type = RAND -- specify number of solutions to generate
$nsol    =  100000;     # full moment tensor
$nsol_fixlam =  100000;     # fixed lambda (includes DC)

# number of entries for flag R
$nR = 0;    # if not specified then it's the full range

# flag for old grid
$oldgrid = 0;   # default is oldgrid=0 (ie do not use old grid)

#----------------------------------------------------------- 

# number of freedom per sample for estimating uncertainty
$nof = 0.01;
$dt = 0.1;

# rms thresholds for discarding bad traces
@thrshd = (10., 10., 10., 10., 10.);

# paramters for running cap for multiple depths
$dep_inc = 0;

# command line input, [] means optional, see () for default value
$usage = 
"
Usage: cap.pl -Mmodel_depth/mag 
    [-A<dep_min/dep_max/dep_inc>] [-C<f1_pnl/f2_pnl/f1_sw/f2_sw>]
    [-D<w1/p1/p2>] [-F<thr>] [-Ggreen] [-Hdt] 
    [-I<nsol> OR -I<Nv/Nw/Nstrike/Ndip/Nrake>]
    [-k1 (old grid setup)] [-L<tau>]
    [-M$model_$dep] [-m$mw OR -m<mw1>/<mw2>/<dmw> ] [-N<n>]
    [-O] [-P[<Yscale[/Xscale_b[/Xscale_s[/k]]]]>] [-Qnof]
    [-R<v0/w0/strike0/dip0/rake0> OR -R<v1/v2/w1/w1/strike1/strike2/dip1/dip2/rake1/rake2>] 
    [-S<s1/s2[/tie]>] [-T<tmax_body/tmax_surf>]
    [-Udirct] [-V<vp/vl/vr>] [-Wi] [-Y<norm>] [-Zstring] event_dirs

    -A  run cap for different depths. (dep_min/dep_max/dep_inc).
    -B  apply same surface wave static shift to all stations. When used static_shift from weight file (last column) is ignored.
    -C  filters for Pnl and surface waves, specified by the corner
        frequencies of the band-pass filter. ($f1_pnl/$f2_pnl/$f1_sw/$f2_sw).
    -D	weight for Pnl (w1) and distance scaling powers for Pnl (p1) and surface waves (p2).
        If p1 is 0 (or not an integer), then Pnl will be turned off.
        Currently there is no option to turn off the surface waves with this command
        or to upweight the surface waves; this must be done via the weight file.
        If p1 or p2 is negative, all traces will be normalized.
        ($weight_of_pnl/$power_of_body/$power_of_surf)
    -E  To up weight rayleigh and surface waves ($weight_of_rayleigh/$weight_of_love)
    -F	include first-motion data in the search. thr is the threshold ($fm_thr).
    	The first motion data are specified in $weight. The polarities
        can be specified using +-1 for P, +-2 for SV, and +-3 for SH after
        the station name, e.g. LHSA/+1/-3 means that P is up and SH is CCW.
        The Green functions need to have take-off angles stored in the SAC
        header.
        threshold should be NEGATIVE if polarities are allowed to conflict expected polarity (as mentioned in weight file).
    -G  Green's function library location ($green).
    -H  dt ($dt).
    -I  specify number of solutions (random mode) OR number of points per parameter.
        RAND: -I<nsol>  e.g. -I10000  --- will generate 10,000 random solutions.
        GRID: -I<Nv>/<Nw>/<Nstrike>/<Ndip>/<Nrake> where Nx = number of poits for parameter x
    -J  Specify cap version to run, eg: `cap.pl -Jlala ...` expects a binary called `caplala`
        This is useful when using multiple cap versions (eg one that writes a binary file)
        If not specified this defaults to `cap`.
    -K  sort output file by distance (=1) or azimuth (=2). default (=0) does nothing.
    -k  Specify k1 to build a lune grid, not a UV grid.
        Use only when LUNE_GRID_INSTEAD_OF_UV = 0 (and recompile cap)
    -L  source duration (estimate from mw, can put a sac file name here).
    -M	specify the model and source depth.
    -m	specify point magnitude: -m<mw0> OR magnitude range: -m<mw1>/<mw2>/<dmw>
    -N  repeat the inversion n times and discard bad traces ($repeat).
    -O  output CAP input (off).
    -P	generate waveform-fit plot with plotting scale. ([-P<Yscale>/<XscaleBody>/<XscaleSurf>](/k))
    	Yscale: amplitude in inch for the first trace of the page ($amplify).
        # scale down: -P1e-5, -P1e-6,... seismograms scaled by this amplitude (smaller the number larger the waveform)
        # scale up:   -P1, -P2,...       seismograms also scaled by amplitude but then enlarged by factor 1, 2, (anything greater than 1)... (larger the number larger the waveform)
        # scale each window:  -P0.5e+0.5 seismograms will be scaled for each component window (all same size) - Normalized scaling
        Xscale: seconds per inch. (body: $spib, surface:$spis).
        append k if one wants to keep those waveforms.
        -p  (small p) For amplitude scaling of surface waves; example: If set to 2 surface waves amplitude will be multipled by twice 
    -Q  number of freedom per sample ($nof)
    -R	For double couple use -R0/0.
        For point solution use -Rv0/w0/strike0/dip0/rake0
        For search range use -Rv1/v2/w1/w2/strike1/strike2/dip1/dip2/rake1/rake2
        Note: This should come after -I flag in the command or it crashes sometimes! 
    -S	max. time shifts in sec for Pnl and surface waves ($max_shft1/$max_shift2) and
        tie between SH shift and SV shift:
        tie=0 		shift SV and SH independently,
        tie=0.5 	force the same shift for SH and SV ($tie).
    -T	max. time window lengths for Pnl and surface waves ($tmax_body/$tmax_surf).
    -U  directivity, specify rupture direction on the fault plane (off).
    -V	apparent velocities for Pnl, Love, and Rayleigh waves (off).
    -W  Integration.
        W0 => Do not integrate. Uses data in its original form
        W1 => (default) Integrate
    -X  weight for normalized polarity misfit [0,1). This is combined with the waveform misfit.
        Recommended value: X0.5 (Bug: Crashes when using -X flag and no polarity is available in the weight file)
        Not recommended: X1. This allows multiple solutions that could fit the observed polarities.
        use -X.99 instead in order to include at least some waveform measure and find a single best solution.

    -Y  specify norm (1 - L1 norm; 2 - L2 norm)
    -Z  specify a different weight file name ($weight).

For details about data preparation and usage, see ``readme_01_preparation.md`` and ``readme_02_usage.md``
";

@ARGV > 0 || die $usage;

$ncom = 5;	# 5 segments to plot

open(INP,">$inp_cmd");
print INP "cap.pl ";
foreach $argnum (0 .. $#ARGV) {
    print INP "@ARGV[$argnum] ";}
close(INP);
system("chmod +x $inp_cmd");

#input options
foreach (grep(/^-/,@ARGV)) {
   $opt = substr($_,1,1);
   @value = split(/\//,substr($_,2));
   if ($opt eq "A") {
     $dep_min = $value[0];
     $dep_max = $value[1];
     $dep_inc = $value[2];
     if ($#value ==2){
       printf STDERR "Running cap for multiple depths: $dep_min to $dep_max at $dep_inc km increment\nWarning: overwriting -Mdepth\n";
     } else {
       $dep_inc=0;
       printf STDERR "Depth run flag -A not specified correctly\nUsing -Mdepth instead\n---------------------\n"; }
   } elsif ($opt eq "B") {
     ($Sstatic_shift) = @value;
   } elsif ($opt eq "C") {
     ($f1_pnl, $f2_pnl, $f1_sw, $f2_sw) = @value;
   } elsif ($opt eq "D") {
     ($weight_of_pnl,$power_of_body,$power_of_surf) = @value;
   } elsif ($opt eq "E") {
     ($weight_of_rayleigh,$weight_of_love) = @value;
   } elsif ($opt eq "F") {
     $fm_thr = $value[0] if $#value >= 0;
   } elsif ($opt eq "G") {
     $green = substr($_,2);
   } elsif ($opt eq "H") {
     $dt = $value[0];
   } elsif ($opt eq "I") {
       # Parameter for number of points or number of solutions
       # Two options.
       #    Length 1 = (RAND) number of solutions (nsol)
       #    Length 5 = (GRID) nv/nw/nk/nh/ns. nsol = nv*nw*nk*nh*ns
       #    
#      $deg = $value[0];
#      $dm = $value[1] if $#value > 0;
#      $dlune = $value[2] if $value[2]
       if ($#value==0) {
           $nI = 1;
           $nsol = $value[0];
           $nv = $nw = $nk = $nh = $ns = $nsol;
       } elsif ($#value==4) {
           $nI = 5;
           ($nv, $nw, $nk, $nh, $ns) = @value;
           $nsol = $nv * $nw * $nk * $nh * $ns;
           if ($nsol <= 0) {
               print STDERR "ERROR. Number of points per paramete should be >0.\n";
               exit(0);
           }
       }
#  elsif ($opt eq "J") {
#    $iso1   = $value[0] if $value[0];
#    $iso2  = $value[1] if $value[1];
#    $clvd1  = $value[2] if $value[2];
#    $clvd2 = $value[3] if $value[3];
#    $fmt_flag="true";     # used later for renaming output figures with fmt
#  }
   } elsif ($opt eq "J") {
     $capver = @value[0];
   } elsif ($opt eq "K") {
     $isort = $value[0];
   } elsif ($opt eq "k") {
     $oldgrid = $value[0];
   } elsif ($opt eq "L") {
     $dura = join('/',@value);
   } elsif ($opt eq "m") {
       # Search parameter for magnitude
       # Two options.
       # 1. Length 1 = fixed magnitude 
       # 2. Length 3 = search over magnitude range 
       # 
       #($md_dep,$mg) = @value;
       #
       if ($#value==0) {
           ($mw1, $mw2, $dm) = ($value[0], $value[0], 0);
           $mg = $value[0];     # IS THIS NEEDED ANYMORE?
       } elsif ($#value==2) {
           ($mw1, $mw2, $dm) = @value;
           # nmw = number of magnitude points when running magnitude search.
           # note this allows to run at grid points dmw finer than 
           # 0.1 (eg 0.01, 0.001...), so run with care
           $nmw = sprintf("%.0f", ($mw2 - $mw1) / $dm +1 );
       }
   } elsif ($opt eq "M") {
       # ($md_dep,$mg) = @value;
       $md_dep = @value[0];
   } elsif ($opt eq "N") {
        $repeat = $value[0];
   } elsif ($opt eq "O") {
     $cmd = "cat";
   } elsif ($opt eq "P") {
     $plot = 1;
     $ampbody_input = $value[0] if $#value >= 0;
     $spib = $value[1] if $value[1] > 0;
     $spis = $value[2] if $value[2] > 0;
     $keep = 1 if $#value > 2;
   } elsif ($opt eq "p") {
     $ampsurf_input = $value[0];
   } elsif ($opt eq "Q") {
     $nof = $value[0];
   } elsif ($opt eq "R") {
       # Flag to set Ranges of parameters
       # Four options: 
       # 1. no flag  = FMT over full range
       # 2. Length 2 = fixed eigenvalue with grid search (v0/w0)
       # 3. Length 5 = fixed moment tensor (v0/w0/strike0/dip0/rake0)
       # 4. Length 10 = subset case (v1/v2/w1/w2/strike1/strike2/dip1/dip2/rake1/rake2)
       if ($#value==1) {
           $nR = 2;
           ($v1, $v2) = @value;
           ($w1, $w2) = @value;
           $nv = $nw = 1;       # at least one lune point needs to exist
       } elsif ($#value==4) {
           $nR = 5;
           ($v1, $w1, $k1, $h1, $s1) = @value;
           ($v2, $w2, $k2, $h2, $s2) = @value;
           $h1 = $h2 = cos($value[3]*$deg2rad);  # read dip from command line and converts on-the-fly to h
           #$h1 = $h2 = $value[3];         # read h from command line
           $nsol = $nv = $nw = $nk = $nh = $ns = 1;
       } elsif ($#value==9) {
           $nR = 10;
           ($v1, $v2, $w1, $w2, $k1, $k2, $dip1, $dip2, $s1, $s2) = @value;
           $h1 = cos($dip1*$deg2rad);
           $h2 = cos($dip2*$deg2rad);
       }
   } elsif ($opt eq "S") {
     ($max_shft1, $max_shft2) = @value;
     $tie = $value[2] if $#value > 1;
   } elsif ($opt eq "T") {
     ($tmax_body, $tmax_surf) = @value;
   } elsif ($opt eq "U") {
     ($rupDir) = @value;
     $pVel = 6.4;
     $sVel = 3.6;
     $rise = 0.4;
     $dirct = "_dir";
   } elsif ($opt eq "V") {
     ($vp, $love, $rayleigh) = @value;
   } elsif ($opt eq "W") {
     $disp = $value[0];
   } elsif ($opt eq "X") {
     $pol_wt = $value[0];
   } elsif ($opt eq "Y") {
     $norm = $value[0];
   } elsif ($opt eq "Z") {
     $weight = $value[0];
   } else {
     printf STDERR $usage;
     exit(0);
   }
}
@event = grep(!/^-/,@ARGV);

#-----------------------------------------------------------
#  Read parameter file instead
if (@ARGV == 1) {
   $parameter_file = $ARGV[0];
   sub_read_parameter_file($parameter_file)
}
#-----------------------------------------------------------

#-----------------------------------------------------------
# prepare additional parameters for CAP input
#-----------------------------------------------------------

# start to output some values
print STDERR "-------------------------------------------------------------\n";
print STDERR "cap.pl: input parameters for CAP:\n";

# convert strike/dip/rake to radians
if ($oldgrid == 0) {
    $k1 = $k1 * $deg2rad;
    $k2 = $k2 * $deg2rad;
    $s1 = $s1 * $deg2rad;
    $s2 = $s2 * $deg2rad;
}

unless ($dura) {
  $dura = int(10**(($mg-5)/2)+0.5);
  $dura = 1 if $dura < 1;
  $dura = 9 if $dura > 9;
}

# Flag to create regular grid as in Alvizuri & Tape (2016) and Silwal & Tape (2016).
# NOTE reproducibility may not be exact since grid spacing uses function gridvec (in the c code).
# Function gridvec does not implement the discretization of the previous version of 
# cap.c which uses rules to account for special grid points.
# Function gridvec also avoids endpoints in all parameters.
if( ($oldgrid == 1) && ($nI == 5)) {
    # Old grid
    print STDERR "\nWarning. Using flag -k1 for the lune grid. cap.c should\n";
    print STDERR "\tbe compiled with flag LUNE_GRID_INSTEAD_OF_UV = 1\n\n";
    # check that K flag works with flag R
    if (($nv == 1) && ($nw == 1) && ($nR == 0)) {
        # default is double couple if Range not specified and it's a single lune point
        ($dv, $dw) = (0, 0);
        ($v1, $v2) = (0, 0);
        ($w1, $w2) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
    } elsif (($nv == 1) && ($nw == 1) && ($nR == 5)) {
        print STDERR "Fixed solution. Input values: $v1 $w1 $k1 $h1 $s1\n";
        # # if Range is set then its a subset
        # # lune points come from user input
        ($dv, $dw) = (0, 0);
    } elsif (($nv == 1) && ($nw == 1) && ($nR > 5)) {
        # If range=5 then it's a point solution.
        # Lune poionts come from user input
        ($dv, $dw) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
    } else {
        # set the full range
        ($v1, $v2) = (-30, 30);
        ($w1, $w2) = (-90, 90);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
    }
    # NOTE values nX are grid spacings, NOT number of points
    # NOTE values dX should be integers (CAP expects integers)
    $dk = sprintf("%.0f", (($k2 - $k1) / $nk));  # strike -- do not include 360
    $dh = sprintf("%.0f", (($h2 - $h1) / $nh));  # dip    -- include 0 at start (though it will be offset later)
    $ds = sprintf("%.0f", (($s2 - $s1) / $ns));  # rake   -- include 0 point

    #print STDERR "Input parameters:\n$dv $dw $dk $dh $ds\n";  # for debugging
    $nsol = $nv * $nw * $nk * $nh * $ns;
} elsif( ($oldgrid == 1) && ($nI == 1)) {
    # random mode
    print STDERR "Warning. Using the old random setup.\n";
    print STDERR "cap.c should be compiled with flag LUNE_GRID_INSTEAD_OF_UV = 1\n";
    if (($nv == 1) && ($nw == 1) && ($nR == 0)) {
        # default is double couple if Range not specified and it's a single lune point
        ($dv, $dw) = (0, 0);
        ($v1, $v2) = (0, 0);
        ($w1, $w2) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
        $nk = $nh = $ns = $nsol;
    } elsif (($nv == 1) && ($nw == 1) && ($nR == 5)) {
        print STDERR "Fixed solution. Input values: $v1 $w1 $k1 $h1 $s1\n";
        # # if Range is set then its a subset
        # # lune points come from user input
        ($dv, $dw) = (0, 0);
    } elsif (($nv == 1) && ($nw == 1) && ($nR > 0)) {
        # if Range is set then its a subset
        # lune points come from user input
        ($dv, $dw) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
        $nk = $nh = $ns = $nsol;
    } else {
        # set the full range
        ($v1, $v2) = (-30, 30);
        ($w1, $w2) = (-90, 90);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
        $nv = $nw = $nk = $nh = $ns = $nsol;
    }
}

# plots for the DC don't have "fmt" in their filenames
if (($v1 == 0) && ($w1 == 0)) {
    $fmt_flag="false";   # double couple
    print STDERR "NOTE computing a double couple solution\n";
} else {
    $fmt_flag="true";
}

#-----------------------------------------------------------
# CHECK THAT USER INPUT MAKE SENSE
#-----------------------------------------------------------

# Set grid type to build (uniform rand, uniform grid) depending on entries in I (previously flag K).
# If using flag k3 then also need to set LUNE_GRID_INSTEAD_OF_UV=1 in cap.c
if ($nI == 1) {
    $grid_type = 2; # RAND
    $grid_type_label="random"
} elsif ($nI == 5) {
    $grid_type = 1; # GRID
    $grid_type_label="grid"
} else{
    print STDERR "STOP. Unable to set search type, check entries in flag R or I\n";
    exit(0);
}

#$grid_type = $type;

# Flag I: set defaults
if (($nI == 5) && ($nv == 1) && ($nw == 1) && ($oldgrid == 0)) {
    # default to the double couple
    ($v1, $v2) = (0, 0);
    ($w1, $w2) = (0, 0);
}

if ( -r $dura ) {	# use a sac file for source time function   
  $dt = 0;
  $riseTime = 1;
} else {
  $riseTime = $rise*$dura;
}

($model, $depth) = split('_', $md_dep);
unless ($depth) {
  $model = ".";
  $depth = 1;
}

if ($dep_inc==0) {
  $dep_min=$depth;
  $dep_max=$depth;
  $dep_inc=1;
}

# convert filter frequencies to periods
$Tf1 = 1/$f1_pnl;
$Tf2 = 1/$f2_pnl;
$Tf3 = 1/$f1_sw;
$Tf4 = 1/$f2_sw;
$filterBand = sprintf("Body:%.2f-%.2f. Surf:%.2f-%.2f",$Tf2,$Tf1,$Tf4,$Tf3);

#-----------------------------------------------------------
# END prepare additional parameters for CAP input
#-----------------------------------------------------------

#-----------------------------------------------------------
# LOOP INVERSIONS BY DEPTH
#-----------------------------------------------------------
for($dep=$dep_min;$dep<=$dep_max;$dep=$dep+$dep_inc) {
    foreach $eve (@event) {

        $md_dep = $model.'_'.$dep;
        next unless -d $eve;
        print STDERR "EVENT ID = $eve | EVENT DEPTH = $dep |  SOURCE DURATION = $dura\n";
        print STDERR "GRID TYPE = $grid_type ($grid_type_label search) | NSOL = $nsol\n";
        print STDERR "-------------------------------------------------------------\n\n";

        #-----------------------------------------------------------
        # copy, sort, cleanup input weight file
        #-----------------------------------------------------------
        # Sort weight file by distance or azimuth
        # Default is to use the weight file as it is
        $input_weight_file = "$eve/$weight";
        $clean_weight_file = "$eve/WEIGHT_CLEAN.dat";
        $station_info_file = "$eve/${eve}_station_list_ALL.dat";

        print STDERR "Sorting and cleaning up weights file: $input_weight_file ...\n";
        open(IN,$input_weight_file) || die "STOP. Weight file not found: $input_weight_file\n";
        @weightlines = <IN>; $nsta_wfile = @weightlines;
        close(IN);

        # THIS SECTION DOES NOT RUN IF K IS NOT SPECIFIED
        # BUG: VR changes whether K is used or not.
        #      To avoid this we might have to run this section by default.
        if ($isort != 0) {
            if ($isort == 1){
                $isortstr = "distance";
            }
            elsif ($isort == 2){
                $isortstr = "azimuth";
            }
            else {
                die "STOP. Sorting can only be by distance (-K1) or azimuth (-K2)\n";
            }

            # Get data from station_list_ALL.dat
            # sta net   lat     lon         dist      azim
            # BKS BK 37.876200 -122.235600 518.062577 279.769407
            # CMB BK 38.034500 -120.386500 360.643229 285.609184
            # MHC BK 37.341600 -121.642600 462.460019 273.168208
            open(IN, $station_info_file) || die "STOP. Station file not found: ${eve}/${eve}_station_list_ALL.dat\n";
            @stnlines = <IN>;
            $nsta_infofile = @stnlines;
            close(IN);

            if ($nsta_wfile > $nsta_infofile) {
                die "STOP. Mismatch NSTA. Not enough information to sort. wfile: $nsta_wfile. station file: $nsta_infofile\n";
            }
            elsif ($nsta_wfile < $nsta_infofile) {
                print STDERR "WARNING. Mismatch NSTA. Check that wfile has the right station listing.\nwfile: $nsta_wfile. station file: $nsta_infofile\n";
            }

            # Read specific columns from the weight file
            print STDERR "Sorting stations by $isortstr ...";
            for ($i = 0; $i < $nsta_wfile; $i++) {
                # NOTE ``split`` does not complain if the weight file doesn't have a column for shift2. 
                # split will still parse the columns correctly.
                ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift,$shift2)=split(" ",@weightlines[$i]);
                ($stnm,$pol) = split("/",$name);
                ($eve1,$net1,$name1,$loc1,$chan1) = split(/\./,$stnm);

                # sort by distance, azimuth, or leave as is
                for ($j = 0; $j < $nsta_infofile; $j++) {
                    # read columns from station_list_ALL.dat
                    ($name2, $net2, $lat2, $lon2, $dist2, $az2) = split(" ",@stnlines[$j]);
                    if (($name1 eq $name2) && ($net1 eq $net2)) {
                        # by distance
                        if ($isort == 1) {
                            $col_to_sort[$i] = $dist2;
                        }
                        # by azimuth
                        elsif ($isort == 2) {
                            $col_to_sort[$i] = $az2;
                        }
                        # THIS CODE IS NEVER REACHED
                        # do nothing
                        else {
                            $col_to_sort[$i] = $i;
                        }
                    }
                }
            }
            print STDERR "done\n";
            # GET INDICES OF SORTED ELEMENTS IN station_list_ALL.dat
            # https://www.perltutorial.org/perl-sort/
            # https://perldoc.perl.org/functions/sort
            # <=>: spaceship operator. sorts numerically.
            # $a, $b are implicitly local, pre defined in perl's sort routine.
            print STDERR "Station distances: @col_to_sort\n";
            @sorted_indices = sort{$col_to_sort[$a] <=> $col_to_sort[$b]} 0 .. $#col_to_sort;

            #-----------------------------------------------------------
            # Remove stations that have no information. 
            # DO THIS REGARDLESS. AT PRESENT THE SCRIPT SEEMS TO TAKE A STATION
            # INTO ACCOUNT WHETHER THE STATION IS USED OR NOT
            #-----------------------------------------------------------
            print STDERR "Removing unused rows ...\n";
            # NO polarity, NO P weight, NO S weight
            open(OUT,'>',$clean_weight_file);
            for ($i = 0; $i < $nsta_wfile; $i++) {
                $ipol = 1;
                # Copy tshifts surf-->love if love empty.
                $ncol_weight = split(" ",@weightlines[$sorted_indices[$i]]);
                if ($ncol_weight == 12) {
                    print STDERR "WARNING. Weight input has 12 columns (OLD). ";
                    ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift) = split(" ",@weightlines[$sorted_indices[$i]]);
                    $shift2 = $shift;
                    print STDERR "Copying col12-->col13 (surf-->love). tshift $shift2 sec (CHECK!)\n";
                } elsif ($ncol_weight == 13) {
                    # **** CHECK HERE ****
                    #   NOTE sorted_indices comes from station_list_ALL.dat, not from the weight files.
                    #   NOTE the following sort is done on the weight files, but the
                    #        sort index comes from a different file which may be sorted
                    #        differently. Does this always work?
                    ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift,$shift2) = split(" ",@weightlines[$sorted_indices[$i]]);
                } else {
                    die "STOP. Weight file has $ncol_weight columns. Expecting 12 or 13.\n";
                }

                ($stnm,$pol) = split("/",$name);

                # no polarity information
                if ($pol eq ''){
                    $ipol = 0;
                    #print STDERR "sta $stnm ipol $ipol\n";
                    #die "Into this loop!";
                }  
                # If pol or weight checks fail then keep the station
                # Skip IF and(pol checks) and(weight checks)
                # 2022-05-03 CHECK this portion doesn't seem to always work
                #print STDERR "\nDEBUG sta $stnm pol $pol polwt $pol_wt pv $pv pr $pr sv $sv sr $sr st $st >> ";
                if (($ipol==0 || $pol_wt==0 || $pol_wt==999) && $pv==0 && $pr==0 && $sv==0 && $sr==0 && $st==0) {
                    print STDERR "DISCARDING: $name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \t $shift2 \n";
                    #die "Into this loop!";
                    next;
                } 
                # save in new weight file
                else {
                    print OUT "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \t $shift2 \n";
                }
            }
            close(OUT);

            ## Clean this weight file by removing stations that have no information (no polarity, no P weight, no S weight)
            #open(OUT,'>',$clean_weight_file);
            #@sorted_indices = sort{$col_to_sort[$a] <=> $col_to_sort[$b]} 0 .. $#col_to_sort;
            #for ($i = 0; $i < $nsta_wfile; $i++){
            #    $ipol = 1;
            #    ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift)=split(" ",@weightlines[$sorted_indices[$i]]);
            #    ($stnm,$pol) = split("/",$name);
            #
            #    if ($pol eq ''){$ipol = 0;}  # no polarity information
            #    if (($ipol==0 || $pol_wt==0 || $pol_wt==999) && $pv==0 && $pr==0 && $sv==0 && $sr==0 && $st==0){
            #        print STDERR "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \n";
            #        next;} # No information available - skip this station
            #    else {     # save in new weight file
            #        print OUT "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \n";
            #    }
            #}
            #close(OUT);
            print STDERR "\nDone sorting/cleaning weight file. Outfile: $clean_weight_file\n";
        }
        else {
            copy $input_weight_file, $clean_weight_file;
        }
        # --------------------------

        open(WEI,$clean_weight_file);
        @wwf=<WEI>;
        close(WEI);
        $ncom = 2 if $wwf[0] =~ / -1 /;

        $cmd = "cap$capver$dirct $eve $md_dep" unless $cmd eq "cat";
        open(SRC, "| $cmd") || die "can not run $cmd\n";
        print SRC "$pVel $sVel $riseTime $dura $rupDir\n",$riseTime if $dirct eq "_dir";
        print SRC "$model $dep\n";          # first input in regular cap run
        print SRC "$tmax_body $tmax_surf $max_shft1 $max_shft2 $repeat $fm_thr $tie $Sstatic_shift\n";
        print SRC "@thrshd\n" if $repeat;   # no value in regular cap run
        print SRC "$vp $love $rayleigh\n";  # vp, vs1, vs2 (in cap.c)
        print SRC "$power_of_body $power_of_surf $nof\n";
        print SRC "$weight_of_pnl $weight_of_rayleigh $weight_of_love\n";
        print SRC "$plot\n";
        print SRC "$disp $pol_wt\n";
        print SRC "$green/$model/\n";
        print SRC "$grid_type\n";
        print SRC "$norm\n";
        print SRC "$dt $dura $riseTime\n";
        print SRC "$f1_pnl $f2_pnl $f1_sw $f2_sw\n";
        print SRC "$mw1 $mw2 $nmw $dm\n";
        print SRC "$v1 $v2 $nv $dv\n";
        print SRC "$w1 $w2 $nw $dw\n";
        print SRC "$k1 $k2 $nk $dk\n";
        print SRC "$h1 $h2 $nh $dh\n";
        print SRC "$s1 $s2 $ns $ds\n";
        print SRC "$nsol\n";
        printf SRC "%d\n",$#wwf + 1;
        print SRC @wwf;
        close(SRC);
        print STDERR ">> cap.pl: end grid search.\n";

        #-----save a copy of inpur command and weight file in the OUTPUT_DIR
        system("cp", $input_weight_file, './OUTPUT_DIR/weight.dat');
        system("cp", $inp_cmd, "./OUTPUT_DIR/${eve}_${model}_${dep}_caprun");
        #system("git log | head -12 > ./OUTPUT_DIR/last_2git_commits.txt");

        plot:
        if ( $plot > 0 && ($? >> 8) == 0 ) {
            #print STDERR "----------------------------------\n";
            #print STDERR ">> cap.pl: plot results ... \n";
            $odir = "./OUTPUT_DIR";
            chdir($odir);
            @dum = split('_', $md_dep);  # split mdl string
            $outfile = sprintf("%s_%s_%03d.out", @event, $model, int($dep));
            open(my $out,'>>',$outfile);
            say $out "INPUT_PAR $md_dep P_win $tmax_body S_win $tmax_surf P $ampbody_input p $ampsurf_input NCOM $ncom spiB $spib spiS $spis $filterBand FMT $fmt_flag";
            print STDERR ">> cap.pl: spis $spis S_win $tmax_surf \n";

            &plot($md_dep, $tmax_body, $tmax_surf, $ampbody_input, $ampsurf_input, $ncom, $spib, $spis, $filterBand, $fmt_flag, @event, $model, $dep, $dura, $riseTime, $pol_wt);
            unlink(<${md_dep}_*.?>) unless $keep;
            chdir("../");
            print STDERR "cap.pl: Inversion and plotting completed.\n";
        } else {
            print STDERR ">> cap.pl: no plots generated.\n";
        }
    }
}
# rename run command and output dir
(my $OUTDIR      = sprintf("fmtout_@{event}_${model}_${dep}_runtime_${runtime1}"));
(my $run_command = sprintf("fmtrun_@{event}_${model}_${dep}_runtime_${runtime2}")); # append the following to supress newline: ``=~ s/\s//g;``
open(INP,">$run_command");
print INP "cap.pl ";
foreach $argnum (0 .. $#ARGV) {
    print INP "@ARGV[$argnum] ";
}
print INP "\n";
close(INP);
system("chmod +x $run_command");
rename("$inp_cmd", $run_command);
#rename("OUTPUT_DIR", "$OUTDIR") || die ( "Error in renaming" );
#print STDERR "cap.pl: Output dir: $OUTDIR\n";

exit(0);
