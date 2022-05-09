# this subroutine plots waveform fits produced by source inversion srct
use List::Util qw[min max];

sub plot {
    print STDERR "\n============================\n";
    print STDERR "cap_plt: Begin plotting results ...\n";
    print STDERR "============================\n";

#  local($mdl, $tmax_body, $tmax_surf $ampbody_input, $num_com, $spis) = @_; # original
  local($mdl, $tmax_body, $tmax_surf, $ampbody_input, $ampsurf_input, $num_com, $spib, $spis, $filterBand, $fmt_flag, $evid, $model, $depth, $dura, $riseTime, $pol_wt) = @_;
  local($nn,$tt,$plt1,$plt2,$plt3,$plt4,$i,$nam,$com1,$com2,$j,$x,$y,@aa,$capout,@name,@aztk);

# set this =1 if you want to plot time windows that have been excluded
  local $keepBad = 0;

# if you want to plot only polarities on the big beachball plot (No azimuth, station info or title)  - default 0
$only_pol = 0;
  
  @trace = ("1/255/255/255","3/0/0/0");       # plot data trace
  @name = ("P V","P R","Surf V"," Surf R","Surf T");

  $filterBand = "Filter periods (seconds): $filterBand";    # 20120719 - report filter bands
  $dura = sprintf("%.2f",$dura);
  $riseTime = sprintf("%.2f",$riseTime);
  $duration = "duration: $dura/$riseTime s";
  
#--------------------------
# plotting GMT file with synthetics and data
  
  # set GMT defaults
  # for all options (such as PAPER_MEDIA): http://gmt.soest.hawaii.edu/gmt/html/man/gmtdefaults.html
  # to check defaults, e.g.: gmtdefaults -L | grep MEASURE_UNIT
  @dum = split('_', $mdl);  # split mdl string
  $ftag=sprintf("%s_%s_%03d", $evid, $model, int($depth));
  $capout_file = sprintf("%s.out", $ftag);

#  $outps = "$mdl.ps";   # original
  #$outps = sprintf("%s_%s_%03d.ps", $evid, $model, int($depth)); # reformatted filename
  #$outps = sprintf("%s_%s_%03d_fmt.ps", $evid, $model, ,int($depth)) if $fmt_flag eq "true";
  $outps = sprintf("${ftag}.ps");
  $outps = sprintf("${ftag}_fmt.ps") if $fmt_flag eq "true";

  # read in the output file results
#  open(FFF,"$mdl.out"); # original
  open(FFF,$capout_file);    # 20130102 calvizuri - new file name
  @capout = <FFF>;
  close(FFF);
  @meca     = split('\s+',shift(@capout));   # Line 1
  @hypo     = split('\s+',shift(@capout));   # Line 2
  @tensor   = split('\s+',$capout[0]);     # Line 3
  @others   = grep(/^#/,@capout);          # Line 4
  @ncomp    = grep(/^#/,@capout);           # Line 5
  @capout   = grep(!/^#/,@capout);             # Remaing Lines
  $nrow = @capout;
  printf STDERR "****************** DEBUG tmax_surf $tmax_surf capout $capout \n";

 # check if there are Input parameters in the last line
  @part = ();
  @last=split(' ',$capout[$nrow-1]);
  if ($last[0] eq 'INPUT_PAR') {
      #$nrow = nrow-1;
      for $ii (0..$nrow-2) {
	   push @part, $capout[$ii];}
       @capout=@part;
       $nrow=$nrow-1;
   }

  # Page size
  $pheight_in = $nrow + 2;  # height of pape

  # positions of seismograms on the page
  # height of each seismogram
  $nn = int($pheight_in);
  $height = $pheight_in - 0.5;
  #($nn,$height) = (12,10.5);   # 10 rows of traces per 10.5 in.
  #print "\n$nn rows of traces per $height in";  
  print "\nseconds per inch = $spis";
  $sepb = 0.2*$spib;    # sec per inch (*1/2) bet body waves
  $seps = 0.2*$spis;    # separation bet surface waves

#  ($tt, $inc) = (2*$tmax_body + 3*$tmax_surf + 4*$sepa, 1);
#  ($tt, $inc) = (2*$tmax_body + 3*$tmax_surf + 2*$seps+2*$sepb, 1);
#  ($tt, $inc) = ($tmax_body + $tmax_surf + $sepb, 4) if $num_com == 2;

  ($twin_body, $inc) = (2*$tmax_body + 2*$sepb, 1);
  ($twin_surf, $inc) = (3*$tmax_surf + 2*$seps, 1);  # 2021-04-02 width for surf windows. stretches wiggles, time axis.
  print STDERR "*** DEBUG twin_surf $twin_surf inc $inc\n";
  $tt = $twin_body + $twin_surf;

  ($tt, $inc) = ($tmax_body + $tmax_surf + $sepb, 4) if $num_com == 2;
#  $width = 0.1*int(10*$tt/$sec_per_inch+0.5);
  $widthb = 0.1*int(10*$twin_body/$spib+0.5);
  $widths = 0.1*int(10*$twin_surf/$spis+0.5);
  $width = $widths+$widthb;
#  @x0 = ($tmax_body+$sepa, $tmax_body+$sepa, $tmax_surf+$sepa, $tmax_surf+$sepa, $tmax_surf);
  @x0 = ($tmax_body+$sepb, $tmax_body+$sepb, $tmax_surf+$seps, $tmax_surf+$seps, $tmax_surf);
  print "\n*** DEBUG x0 @x0 *** \n";

  $pwidth_in = $width +1.5 ;  # width of paper    # orig 8.5
  print "\n$nrow rows to plot";
  print "\npaper is $pwidth_in inches wide and $pheight_in inches tall";
  print STDERR ">>> gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch\n";
  system("gmt gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");

  # horizontal offset (why is it needed?)
  #$xoffset="3.0";
  $xoffset=$widthb;

  #===========================================================================================
  # 2022-05-09. THIS SECTION DEALS WITH SCALING AMPLITUDES FOR BODY AND SURF WAVEFORMS.
  # This was a little involved because originally there was a feature (bug?) that allowed to plot all wiggles with the same amplitude.
  # This scaling was used for the Uturuncu paper and needed some fiddling.
  # Then I adapted pssac through various GMT revisions until GMT6 
  # TODO: replace all this with obspy.
  #
  # 2022-05-09 disable the following sections. save for reference Uturuncu paper.
  # KEY: set amplitude scaling for seismograms
  # if ($ampbody_input > 0.)   {$stam = "$ampbody_input/-1";}                                      # 2022-05-09 disabled. not used anywhere else in the code.
  # else                       {$stam = -$ampbody_input;} # original line (with pssac, not pssac2) # 2022-05-09 disabled. not used anywhere else in the code.
  # if ($ampbody_input == 0x0) {$amp = $ampbody_input;}                 # Uturuncu paper to scale all wiggles equally.
  # else                       {$amp = $ampbody_input/$ampsurf_input;}  # Uturuncu paper to scale all wiggles equally.
  # $ampsurf_flag = "$amp/0.";
  # $ampbody_flag = "$ampbody_input/0."; # overwrite for absolute (to match default plotting) # 2022-05-09 disabled. Uturuncu paper to scale all wiggles equally.
  # print "\namplitude scaling ampbody_input = $ampbody_input";
  # #print "\npssac2 amplitude scaling stam = $stam\n";    # 2022-05-09 disabled. not used anywhere else in the code.
  # print "\n*** DEBUG amplitude scaling ampsurf_flag $ampsurf_flag, ampbody_flag $ampbody_flag, ampbody_input $ampbody_input, amp $amp ***\n";

  # print "Max amplitude body waves $ampmax_body \n";
  # print "Max amplitude surf waves $ampmax_surf \n";
##---------------------------------------
# Three options for plotting (and scaling) the waveforms using -P flag (body waves) and -p flag (surface wave)
# Default for both -P and -p flag is 1 (i.e. option 1 in the following comments and using the scaling_factor=1)
# 1. Normalize by maximum body and surface amplitude separatetly, then apply a scaling factor
# 2. Normalized plotting -- data and synthetics have same maximum amplitude for all waveforms
# 3. The default plotting -- scale waveforms by given amplitude
  #----------------------------------------------------------
  # UPDATE 2021-03-29 
  #
  # pssac2 amplitude scaling (GMT4/5?)
  # -M vertical scaling in sacfile_unit/MEASURE_UNIT = size<required> 
  #           size: each trace will normalized to size (in MEASURE_UNIT)
  #               scale =  (1/size) * [data(max) - data(min)]
  #           size/alpha: plot absolute amplitude multiplied by (1/size)*r^alpha
  #               where r is the distance range in km across surface
  #               specifying alpha = 0.0 will give absolute amplitudes
  #               scale = (1/size) * r^alpha
  #           size/s: plot absolute amplitude multiplied by (1/size)*sqrt(sin(gcarc))
  #               where gcarc is the distance in degrees.
  #               scale = (1/size) * sqrt(sin(gcarc))
  #               
  # PSSAC version GMT 6.1.0
  # Flag -M<size>/<alpha>
  # -M Vertical scaling
  #    <size>: each trace will scaled to <size>. The default unit is PROJ_LENGTH_UNIT.
  #       The scale factor is defined as yscale = size*(north-south)/(depmax-depmin)/map_height 
  #    <size>/<alpha>: 
  #       <alpha> < 0, use the same scaling factor for all traces. The scaling factor will scale the first trace to <size>[<u>].
  #       <alpha> = 0, multiply all traces by <size>. No unit is allowed.  (nb 2021-03-31 "specifying alpha = 0.0 will give absolute amplitudes", see above)
  #
  # PSSAC version GMT  6.3
  #   -M<size>/<alpha>
  #    Vertical scaling, with each trace will scaled to <size>. 
  #    The default unit is PROJ_LENGTH_UNIT. 
  #    The scale factor is defined as yscale = size*(north-south)/(depmax-depmin)/map_height. 
  #    Specify <alpha>:
  #    • <alpha> < 0, use the same scaling factor for all traces. The scaling factor will scale the first trace to <size>[<u>].
  #    • <alpha> = 0, multiply all traces by <size>. No unit is allowed.
  #    • <alpha> > 0, multiply all traces by size*r^alpha, r is the distance range in km.
  # 
  #     NB 2022-05-03
  #         want <size> to anything reasonable (1? = no scaling)
  #         want <alpha> negative to scale all by the same amount. unless want different plotting.
  #----------------------------------------------------------

  #      ##OPT 1. scale by the maximum body wave amplitude ($ampmax_body) and then scale by $ampbody_input factor (-P flag)
  #      ##OPT 2. Normalized plotting (using the pssac2 bug) -- data and synthetics have same maximum amplitude
  #      ##OPT 3. default plotting (FUTURE: find a better way to differentiate b/w exponents and rational number) -- scale by given amplitude $ampbody_input  (-P flag)
  #      ## -P flag. body waves
  #      #if    ($ampbody_input>=0.1) {$ampscale_body = $ampmax_body/$ampbody_input;}
  #      #elsif ($ampbody_input==0)   {$ampscale_body = "0.5e+0.5";}
  #      #else                    {$ampscale_body = $ampbody_input;}
  #      ## -p flag. surface waves
  #      #if    ($ampsurf_input>=0.1){$ampscale_surf = $ampmax_surf/$ampsurf_input;}
  #      #elsif ($ampsurf_input==0)  {$ampscale_surf = "0.5e+0.5";}
  #      #else                   {$ampscale_surf = $ampsurf_input;}
  #      #$ampbody_flag = "$ampscale_body/0.";        # set the parameters
  #      #$ampsurf_flag = "$ampscale_surf/0.";
  #      ## 2021-03-29 update for pssac gmt6
  #      # scale BODY_amp. flag -P. -M<size>/<alpha>
  #      # scale SURF_amp. flag -p. -M<size>/<alpha>
  #      # <size>: each trace will scaled to <size>. The default unit is PROJ_LENGTH_UNIT.
  #      # alpha < 0 use the same scaling factor for all traces. The scaling factor will scale the first trace to <size>[<u>]
  #      # alpha = 0 multiply all traces by <size>. No unit is allowed.  (nb 2021-03-31 from above: specifying alpha = 0.0 will give absolute amplitudes)
  #      # alpha > 0 multiply all traces by size*r^alpha, r is the distance range in km
  #      $ampscale_body = "1/-1";
  #      #$ampscale_body = "3000.0/0";
  #      #$ampscale_body = "30/1";
  #      $ampscale_body = "0.1/-1"; # SIMPLEST FLAG. USE. 2021-05-20. Nepal event
  #      $ampscale_body = "0.2/-1"; # SIMPLEST FLAG. USE. 2021-05-20. Nepal event
  #      #$ampscale_body = "0.5/-1"; # SIMPLEST FLAG. USE. 2021-05-20. Nepal event
  #      #$ampscale_body = "0.9/-1"; # SIMPLEST FLAG. USE. 2021-05-20. Nepal event
  #      $ampscale_surf = "1/-1";   
  #      $ampscale_surf = "0.5/-1"; 
  #      #$ampscale_surf = "50/0.0";
  #      #$ampscale_surf = "1/1";   
  #      #
  ## 2022-05-03 TEST AGAIN GMT 6.3. want: <anything/negative>
  $ampscale_body = "1/-1";
  $ampscale_surf = "1/-1";

  print "pssac norm BODY -M<size/alpha>: $ampscale_body \n";
  print "pssac norm SURF -M<size/alpha>: $ampscale_surf \n";
  print "\n*** DEBUG amplitude scaling ampsurf_flag $ampsurf_flag, ampbody_flag $ampbody_flag, ampbody_input $ampbody_input, amp $amp ***\n";
  # 2022-05-09 END REVISED SECTION THAT DEALS WITH SCALING SEISMOGRAM AMPLITUDES
  #===========================================================================================

  # (1) plot cut seismograms with scaled amplitudes (first command: no -O appears)
  $tscale_x = 0.55;
  $tscale_y = $pheight_in - 2.0;
  #$plt1b = "| pssac2 -JX${widthb}i/${height}i -L${spib} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$twin_body/0/$nn -Y0.2i -Ent-2 -M$ampbody_flag -K -P >> $outps";
  #$plt1s = "| pssac2 -JX${widths}i/${height}i -L${spis} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$twin_surf/0/$nn -X${xoffset}i -Ent-2 -M$ampsurf_flag -O -K -P >> $outps";
  #$plt1b = "| pssac -JX${widthb}i/${height}i -S${spib} -M$ampbody_flag -R0/$twin_body/0/$nn               -Y0.2i      -K -P -V >> $outps";
  #$plt1s = "| pssac -JX${widths}i/${height}i -S${spis} -M$ampsurf_flag -R0/$twin_surf/0/$nn -X${xoffset}i          -O -K -P -V >> $outps";
  $plt1b = "| gmt pssac -JX${widthb}i/${height}i -M$ampscale_body -R0/$twin_body/0/$nn               -Y0.2i       -K -P >> $outps";
  $plt1s = "| gmt pssac -JX${widths}i/${height}i -M$ampscale_surf -R0/$twin_surf/0/$nn -X${xoffset}i           -O -K -P >> $outps";

  # (2) plot text labels
  $plt2_stn_info  = "| gmt pstext -JX -R -O -K -N -X-${xoffset}i >> $outps";
  $plt2_wf_info_b = "| gmt pstext -JX${widthb}i/${height}i -R0/$twin_body/0/$nn -O -K -N >> $outps";
  $plt2_wf_info_s = "| gmt pstext -JX${widths}i/${height}i -R0/$twin_surf/0/$nn -X${xoffset}i -O -K -N >> $outps";

  # (3) plot beachballs (solution, followed by possible local minima)
  $ballcolor = "150";
#  $dY = ${pheight_in} - 1.8;    # original
  $dY = ${pheight_in} - 1.6;
  $dX = -0.7-$xoffset;
#  $plt3 = "| psmeca -JX5i/1i -R-1/9/-1/1 -Sa5i -G$ballcolor -Y${dY}i -X-0.7i -O -K >> $outps";
  $plt3 = "| gmt psmeca -JX5i/1i -R-1/9/-1/1 -Sa5i -G$ballcolor -Y${dY}i -X${dX}i -O -K >> $outps";
  $plt3 = "| gmt psmeca -JX5i/1i -R-1/9/-1/1 -Sm8i -G$ballcolor -Y${dY}i -X${dX}i -O -K >> $outps" if $tensor[1] eq "tensor";

  # (4) plot markers on beachball
  # note: -JPa is a basemap for polar coordinates, clockwise from north

  # azimuths
  $plt4b = "| gmt psxy -JPa1i -R0/360/0/1 -Sc0.02i -N -W0.5p,0/0/0 -G255 -O -K >> $outps";

  # supplemental: upper hemisphere piercing points on beachballs (o)
  #$plt4a = "| psxy -JPa1i -R0/360/0/1 -Sc0.08i -N -W0.5p,255/0/0 -O -K >> $outps";

  # default: lower hemisphere piercing points on beachballs (x) (last command: no -K appears)
  $plt4 = "| gmt psxy -JPa1i -R0/360/0/1 -Sx0.10i -N -W0.5p,255/0/0 -G255 -O -K >> $outps";

#  $plt1=$plt2=$plt3="|cat";	# output GMT commands to command window for testing

  # (2.5) plot header information
  # 2021-04-02 calvizuri -- update
  # 
  # pstext [core] 6.1.0 
  #
  # Reads (x,y[,fontinfo,angle,justify],text) from <table> [or stdin].
  # OR (with -M) one or more text paragraphs with formatting info in the segment header.
  # Built-in escape sequences:
  #    @~ toggles between current font and Symbol font.
  #    @%<no>% switches to font number <no>; @%% resets font.
  #    @:<size>: switches font size; @:: resets font size.
  #    @;<color>; switches font color; @;; resets font color.
  #    @+ toggles between normal and superscript mode.
  #    @- toggles between normal and subscript mode.
  #    @# toggles between normal and Small Caps mode.
  #    @_ toggles between normal and underlined text.
  #    @!<char1><char2> makes one composite character.
  #    @. prints the degree symbol.
  #    @@ prints the @ sign itself.
  #    Use @a|c|e|in|o|s|u|A|C|E|N|O|U for accented European characters.
  # (See module documentation for more information).

  $dX = 0.8;
  $dY = 0.3;
  #$plt4_5 = "| pstext -J -R -Y${dY}i -X${dX}i -O -N >> $outps";
  $plt4_5 = "| gmt pstext -J -R -Y${dY}i -X${dX}i -O -N -F+jl+f12p,Helvetica,black >> $outps";

#--------------------------

#  $outps2 = "${mdl}_beach.ps"; # original
  #$outps2 = sprintf("%s_%s_%03d_beach.ps", $evid, $model, int(depth));   # 20130102 calvizuri - revised filename
  #$outps2 = sprintf("%s_%s_%03d_beach_fmt.ps", $evid, $model, int($depth)) if $fmt_flag eq "true";
  $outps2 = sprintf("${ftag}_beach.ps");
  #$outps2 = sprintf("${ftag}_beach_fmt.ps") if $fmt_flag eq "true";
  $outps2 = sprintf("${ftag}_fmt_beach.ps") if $fmt_flag eq "true";

  $fac = 6.5;
  $fac2 = 8.2*$fac;   # original: 5*$fac
  $JP = "-JPa${fac}i";

  # plot beachball
# $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sm${fac2}i -X1i -Y2i -K -P >> $outps2";
  $xplt3 = "| gmt psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sa${fac}i -X1i -Y2i -K -P >> $outps2";
  $xplt3 = "| gmt psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sm${fac2}i -X1i -Y2i -K -P >> $outps2" if $tensor[1] eq "tensor";

  # plot markers on beachball
  # note: -JPa is a basemap for polar coordinates, clockwise from north

  # azimuths
  $xplt4b = "| gmt psxy $JP -R0/360/0/1 -Sc0.02i -N -W0.5p,0/0/0 -G255 -O -K >> $outps2";

  # supplemental: upper hemisphere piercing points on beachballs (o)
  $xplt4a = "| gmt psxy $JP -R0/360/0/1 -Sc0.08i -N -W0.5p,255/0/0 -O -K >> $outps2";

  # default: lower hemisphere piercing points on beachballs (x)
  $xplt4  = "| gmt psxy $JP -R0/360/0/1 -Sx0.10i -N -W0.5p,255/0/0 -G255 -O -K >> $outps2";
  $xplt4c = "| gmt psxy $JP -R0/360/0/1 -St0.30i -N -W1p,0/255/0 -G255 -O -K >> $outps2";  # up polarity (green) - triangle
  $xplt4d = "| gmt psxy $JP -R0/360/0/1 -Si0.30i -N -W1p,0/0/255 -G255 -O -K >> $outps2";  # down polarity (blue) - triangle
  $xplt4e = "| gmt psxy $JP -R0/360/0/1 -St0.30i -N -W1p,255/0/0 -G255 -O -K >> $outps2";  # non-matching polarity red) - triangle
  $xplt4f = "| gmt psxy $JP -R0/360/0/1 -Si0.30i -N -W1p,255/0/0 -G255 -O -K >> $outps2";  # non-matching polarity (red) - triangle

  # plot text labels
  $xplt5a = "| gmt pstext $JP -R0/360/0/1 -N -O -K >> $outps2";
  $xplt5b = "| gmt pstext $JP -R0/360/0/1 -N -O -K >> $outps2";
  $xplt5c = "| gmt pstext $JP -R0/360/0/1 -N -O -K >> $outps2";

  # title (LAST COMMAND: no -K appears)
  $xplt6 = "| gmt pstext -JX -R -N -O -Xa0 -Ya7.5 >> $outps2";

#--------------------------
# FIGURE 1: waveform fits with moment tensor
  print STDERR ">>> cap_plt.pl: plot waveform fits and best mechanism ...\n";
  # 
  print "\n============================\n";
  print "INVERSION RESULTS\n";
  print "============================\n";
  print "\nmeca:  \t@meca";
  print "\ntensor:\t@tensor";
  print "\nothers:\t@others";
  print "\nCAP outfile:\n";
  print "---------------------------------\n";
  print "@capout";
  print "---------------------------------\n";

# get strike dip and rake
  $stk = @meca[0];
  $dip = @meca[1];
  $rak = @meca[2];

  # compute piercing points for beachballs
  $ampmax_body=0; # maximum aplitude for pssac plotting (-P flag) - body
  $ampmax_surf=0; # maximum aplitude for pssac plotting (-P flag) - surface
  $i = 0; $j = 0; $i2 = 0; $j2 = 0;
  $pi = 3.14159265358979323846;
  @tklh=(); @tkuh=(); @staz=(); @az=(); @tklh_useweights=(); @staz_useweights=(); @tkuh_useweights=();
  foreach (@capout) {
      #     1                            2              3     4  5   6      7        8        9   10   11  12    13   14       15       16  17       18    19   20       21       22  23
      #     |                            |              |     |  |    |     |        |        |   |    |   |     |    |        |        |   |        |     |    |        |        |   |
      # 20210314141526689.II.BORG.00.BH  88.1/-0.00 0   0.00  0  0.00 -0.00 3.03e-05 3.04e-05 0   0.00 83  0.00  0.24 1.93e-05 1.51e-05 1   3.83 99  1.90  0.06 2.28e-03 2.14e-03 1   5.46 98  1.90 -0.09 2.08e-03 2.28e-03 1   9.48 99  3.05  0.23 3.67e-03 2.92e-03 0  -0.35
    @aa = split;
    if ($aa[7] >$ampmax_body && $aa[2] !=0){$ampmax_body=$aa[7] ;} # maxamp body vertical. pssac plotting -P flag
    if ($aa[14]>$ampmax_body && $aa[9] !=0){$ampmax_body=$aa[14];} # maxamp body radial
    if ($aa[21]>$ampmax_surf && $aa[16]!=0){$ampmax_surf=$aa[21];} # maxamp surf vertical
    if ($aa[28]>$ampmax_surf && $aa[23]!=0){$ampmax_surf=$aa[28];} # maxamp surf radial
    if ($aa[35]>$ampmax_surf && $aa[30]!=0){$ampmax_surf=$aa[35];} # maxamp surf transverse
    $ifmp[$i] = $aa[37];  # first-motion polarity (input - data)
    $ifmpt[$i] = $aa[38];  # first-motion polarity (theoretical)
    $stnm = $aa[0];                              # station name
    #next if $aa[2] == 0;                        # skip if no body waves
    $x = `saclst az user1 f ${mdl}_$aa[0].0`;    # get the azimuth and P take-off angle
    @dd = @aa;
    @aa = split(' ', $x);       # outputs something like this: wes_1_HOYA.LL.TPH..LH.0 323.513 90.72
    @aa_pre = @aa;
    #print "\n--> saclst az user1 f ${mdl}_$aa[0].0";

    # compute polar coordinates azimuth and radius
    # NOTE this part outputs all azimuths in the weight file, even if the
    # station was not used in the inversion. Unless the input weigh files are
    # pre-sorted and clean.
    # WARNING if the weight files are not clean this line may cause a mismatch
    # between STNAME and AZIM 
    # CHECK
    $az[$i] = $aa[1];

    $azvec[$i] = sprintf("%s\n",$aa[1]);
    $staz[$i] = sprintf("%s %f %s\n",$aa[1],1.1,$stnm);      # station azimuth
    if ($aa[2]>90.) {                                        # upper hemisphere
       $rad = sqrt(2.)*cos($aa[2]*$pi/360);
       $tkuh[$j] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);
       $j++;
       # project piercing point to lower hemisphere
       $aa[1] += 180;
       $aa[2]=180-$aa[2];
    }
    $rad = sqrt(2.)*sin($aa[2]*$pi/360);
    $tklh[$i] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);        # lower hemisphere
    if (($dd[37] != 0  && $pol_wt != 0) || $dd[2]!=0 || $dd[9]!=0 || $dd[16]!=0 || $dd[23]!=0 || $dd[30]!=0 || $keepBad!=0){
	$tklh_useweights[$i2] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);
	$staz_useweights[$i2] = sprintf("%s %f %s\n",$aa_pre[1],1.1,$stnm);
	if ($aa_pre[2]>90.) {
	    $tkuh_useweights[$j2] = sprintf("%s %f %s\n",$aa_pre[1],$rad,$stnm);
	    $j2++;
	}
	$i2++;
    }
    $i++;
  }
#--------------------------compute pssac plotting info (scaling factor ampmax_body)
  print "============================\n";
  print "Plotting waveforms + best mechanism ...\n";
#---------------------------------------
  # 20151025 cralvizuri - uncomment this command to normalize surf waves
  #                       This is for figures in Uturuncu FMT paper
  #$ampsurf_flag = $ampbody_flag;

  # 2021-04-02 calvizuri -- update for pssac GMT6
  # pssac2 GMT5
  # [-Sshift]
  # -S <seconds> shift the trace by seconds
  # -L <seconds per MEASURE_UNIT> while poltting on maps <required for maps>\
  #         If your seismograms look choppy and pixelated
  #         Check the value of DOTS_PR_INCH in gmtdefaults
  #         and increase the value using gmtset
  #
  # pssac GMT6 
  # -T Time alignment. 
  # [-T[+t<tmark>][+r<reduce_vel>][+s<shift>]]
  #    +t<tmark> align all trace along time mark. Choose <tmark> from -5(b), -4(e), -3(o), -2(a), 0-9(t0-t9).
  #    +r<reduce_vel> reduce velocity in km/s.
  #    +s<shift> shift all traces by <shift> seconds.
  #
  # 2022-05-03 
  # pssac GMT 6.3 
  # -R<west>/<east>/<south>/<north>[+r]
  #
  #$plt1b = "| pssac2 -JX${widthb}i/${height}i -L${spib} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$twin_body/0/$nn -Y0.2i -Ent-2 -M$ampbody_flag -K -P >> $outps";
  #$plt1s = "| pssac2 -JX${widths}i/${height}i -L${spis} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$twin_surf/0/$nn -X${xoffset}i -Ent-2 -M$ampsurf_flag -O -K -P >> $outps";
  #$plt1b = "| pssac -JX${widthb}i/${height}i -S${spib} -M$ampscale_body -R0/$twin_body/0/$nn               -Y0.2i       -K -P >> $outps";
  #$plt1s = "| pssac -JX${widths}i/${height}i -S${spis} -M$ampscale_surf -R0/$twin_surf/0/$nn -X${xoffset}i           -O -K -P >> $outps";
  $plt1b = "| gmt pssac -JX${widthb}i/${height}i -M$ampscale_body -R0/$twin_body/0/$nn               -Y0.2i       -K -P >> $outps";
  $plt1s = "| gmt pssac -JX${widths}i/${height}i -M$ampscale_surf -R0/$twin_surf/0/$nn -X${xoffset}i           -O -K -P >> $outps";
  #$plt1s = "| pssac -JX${widths}i/${height}i -M$ampscale_surf -R0/$twin_surf/0/$nn -X${xoffset}i    -V -Vc -O -K -P >> $outps";
  #print "*** DEBUG plt1s $plt1s\n";

  # remove the file if it exists
  unlink($outps) if -e $outps;
  unlink($outps2) if -e $outps2;

  # save a copy for the second file
  @capout0=@capout;

  while (@capout) {

#    # plot waveforms
#    open(PLT, $plt1);
#    $i = 0;
#    @capout_splice = splice(@capout,0,$nn-2);
#    foreach (@capout_splice) {
#      @aa = split;
#      $nam = "${mdl}_$aa[0].";
#      $x=0;
#      for($j=0;$j<5;$j+=$inc) {
#        $com1=8-2*$j; $com2=$com1+1;
#	if ($aa[4*$j+2]>0) {
#	   printf PLT "%s %f %f 5/0/0/0\n",$nam.$com1,$x,$nn-$i-2;
#           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
#	} elsif ($keepBad) {
#	   printf PLT "%s %f %f 2/0/255/0\n",$nam.$com1,$x,$nn-$i-2;
#           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
#	}
#        $x = $x + $x0[$j];
#      }
#      $i++;
#    }
#    close(PLT);

    # plot waveforms body waves

## PSSAC version GMT 6.1.0
## <saclist> is an ASCII file (or stdin) which contains the name of SAC files to plot and controlling parameters.
##    Each record has 1, 3 or 4 items:  <filename> [<X> <Y> [<pen>]]. 
##    <filename> is the name of SAC file to plot.
##    <X> and <Y> are the position of seismograms to plot on a map.
##       On linear plots, the default <X> is the begin time of SAC file, which will be adjusted if -T option is used, 
##       the default <Y> is determined by -E option.
##       On geographic plots, the default <X> and <Y> are station longitude and latitude specified in SAC header.
##       The <X> and <Y> given here will override the position determined by command line options.
##    If <pen> is given, it will override the pen from -W option for current SAC file only.
##
## -W Set pen attributes [Default pen is default,black]:
##    <pen> is a comma-separated list of three optional items in the order:
##        <width>[c|i|p], <color>, and <style>[c|i|p].
##    <width> >= 0.0 sets pen width (default units are points); alternatively a pen
##              name: Choose among faint, default, or [thin|thick|fat][er|est], or obese.
##    <color> = (1) <gray> or <red>/<green>/<blue>, all in range 0-255;
##              (2) #rrggbb, all in the range 0-255 using hexadecimal numbers;
##              (3) <c>/<m>/<y>/<k> in 0-100% range;
##              (4) <hue>-<sat>-<val> in ranges 0-360, 0-1, 0-1;
##              (5) any valid color name.
##    <style> = (1) pattern of dashes (-) and dots (.), scaled by <width>;
##              (2) "dashed", "dotted", "dashdot", "dotdash", or "solid";
##              (3) <pattern>[:<offset>]; <pattern> holds lengths (default unit points)
##                  of any number of lines and gaps separated by underscores.
##                 The optional <offset> shifts elements from start of the line [0].
## 
    open(PLT, $plt1b);
    $i = 0;
    @capout_splice = splice(@capout,0,$nn-2);
    foreach (@capout_splice) {   # go over each line in .out file
        @aa = split;
	if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
        $nam = "${mdl}_$aa[0].";
        $x=0;
        for($j=0;$j<2;$j+=$inc) {
            $com1=8-2*$j; $com2=$com1+1;    # seismogram extensions (.0, .1, .2...)
            if ($aa[7*$j+2]>0) {
# <saclist> contains SAC files + plotting parameters. See above for formatting instructions.
#           Each record has 1, 3 or 4 items:  <filename> [<X> <Y> [<pen>]].
                #printf "(j=$j) x=$x\t"; # debug
                #printf PLT "%s %f %f 5/0/0/0\n",  $nam.$com1,$x+0,$nn-$i-2;     # data (black)
                #printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x+0,$nn-$i-2;   # synthetic (red)
                printf PLT "%s %f %f 0.8,black\n", $nam.$com1,$x+0,$nn-$i-2;   # data (black)
                printf PLT "%s %f %f 0.8,red\n",   $nam.$com2,$x+0,$nn-$i-2;   # synthetic (red)
            } elsif ($keepBad) {
                printf PLT "%s %f %f 0.8,green\n",$nam.$com1,$x+0,$nn-$i-2;  # bad data (green)
                printf PLT "%s %f %f 0.8,red\n",  $nam.$com2,$x+0,$nn-$i-2;  # synthetic (red)
            }
            $x = $x + $x0[$j];
        }
#                printf "\n"; # debug
        $i++;
    }
    }
    close(PLT);

    #-----------------------------------------------------------
    print "Plotting station data and labels ... \n";
    #-----------------------------------------------------------
    open(PLT, $plt1s);
    $i = 0;
    foreach (@capout_splice) {
        #     1                            2              3     4  5   6      7        8        9   10   11  12    13   14       15       16  17       18    19   20       21       22  23
        #     |                            |              |     |  |    |     |        |        |   |    |   |     |    |        |        |   |        |     |    |        |        |   |
        # 20210314141526689.II.BORG.00.BH  88.1/-0.00 0   0.00  0  0.00 -0.00 3.03e-05 3.04e-05 0   0.00 83  0.00  0.24 1.93e-05 1.51e-05 1   3.83 99  1.90  0.06 2.28e-03 2.14e-03 1   5.46 98  1.90 -0.09 2.08e-03 2.28e-03 1   9.48 99  3.05  0.23 3.67e-03 2.92e-03 0  -0.35
        @aa = split;
	if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
        $nam = "${mdl}_$aa[0].";
        #$x=$x0[1];
        $x=0;   # 0 aligns with text
        for($j=2;$j<5;$j+=$inc) {
            #printf "*** DEBUG j $j x $x\n"; 
            $com1=8-2*$j; $com2=$com1+1;
            if ($aa[7*$j+2]>0) {
                #printf STDERR "cap_plt.pl. DEBUG x %f j %d x0[j] %f TEST %f\n", $x, $j, $x0[$j], $x0[$j]+ $aa[7*$j+5] ;
                # <saclist> contains SAC files + plotting parameters. See above for formatting instructions. 
                #           Each record has 1, 3 or 4 items:  <filename> [<X> <Y> [<pen>]].
                # NOTE the wiggles already are pre-aligned here
                #printf PLT "%s %f %f 5/0/0/0  \n",$nam.$com1,$x,$nn-$i-2;   # data (black)
                #printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;   # synthetic (red)
                printf PLT "%s %f %f 1.0,black\n", $nam.$com1,$x-$aa[7*$j+5],$nn-$i-2;  # data (black)
                #printf PLT "%s %f %f 0.5,red\n",   $nam.$com2,$x,$nn-$i-2;  # synthetic (red)              # 2022-05-04 ORIGINAL
                printf PLT "%s %f %f 0.9,red\n",   $nam.$com2,$x-$aa[7*$j+5],$nn-$i-2;  # synthetic (red)   # 2022-05-04 UPDATE: Include `aa` shift in the synthetics. Why was this was not done originally (GMT 4.5.15-UAF)?
                #printf STDERR "%s OBS %f %f 0.8,black\n", $nam.$com1,$x,$nn-$i-2;  # data (black)
                #printf STDERR "%s SYN %f %f 0.5,red\n",   $nam.$com2,$x,$nn-$i-2;  # synthetic (red)
                #printf stderr "ind 7*$j+2 aa $aa[7*$j+2] $x0[$j]\n";
            } elsif ($keepBad) {
                printf PLT "%s %f %f 0.8,green\n",$nam.$com1,$x,$nn-$i-2;   # bad data (green)
                printf PLT "%s %f %f 0.8,red\n",  $nam.$com2,$x,$nn-$i-2;   # synthetic (red)
            }
            $x = $x + $x0[$j];
            #$x = $x + $x0[$j] + $aa[7*$j+4];
            #$x = $aa[7*$j+2];
            #printf STDERR "cap_plt.pl. DEBUG x %f j %d x0[j] %f aa_ind %f\n", $x, $j, $x0[$j], 7*$j+2, $aa[7*$j+2];
        }
#                printf "\n"; # debug
        $i++;
    }
    }
    close(PLT);

    
    # text labels
#    open(PLT, $plt2);
#    $y = $nn-2;
#    $i=0;
#    foreach (@capout_splice) {
#      @aa = split;
#      $x = 0;
#      printf PLT "%f %f 10 0 0 1 $aa[0]\n",$x-0.8*$spis,$y;            # station label
#      printf PLT "%f %f 10 0 0 1 $aa[1]\n",$x-0.7*$spis,$y-0.2;        # distance_km/overal time shift
#      printf PLT "%f %f 10 0 0 1 %.1f\n",$x-0.7*$spis,$y-0.4,$az[$i];  # azimuth (see az above)
#      $i=$i+1;
#      for($j=0;$j<5;$j+=$inc) {
#	if ($aa[4*$j+2]>0 || $keepBad) {
#
#          printf "(j=$j) x=$x \t ";
#          printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;  # time shift each wave
#          printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;  # correl value
#	}
#        $x = $x + $x0[$j];
#      }
#      $y--;
#    }

    # plot station info
    
    open(PLT, $plt2_stn_info);
    $y = $nn-2;
    $i=0;
    foreach (@capout_splice) {
# Ruler for reading CAP output
#    0          1       2     3   4   5      6      7       8     9    10  11   12    13      14       15   16   17  18   19    20    21        22    23   24  25   26    27    28     29 
#    |          |       |     |   |   |      |      |       |     |     |   |   |     |        |       |    |     |   |    |    |      |        |     |     |   |    |     |    |       |
# PLMK_XP    11.2/0.14  1   0.67 95 -0.08  0.64 8.19e-07 4.32e-07 1   0.79 80 -0.08 -0.09 8.05e-07 8.79e-07 1   3.48 79  1.89  0.84 9.47e-07 4.10e-07 1   4.47 75  1.89  1.24 1.03e-06 2.98e-07 1   3.68 80  0.23  1.61 7.56e-07 1.51e-07  1   0.45
#                                                                                                                                                                                               |     |   |   |      |    |         |      |    |
#                                                                                                                                                                                               30    31 32   33    34    35        36     37   38

        # variables from cap output
        @aa = split;
        @ab = split('/',$aa[1]);
        $dist_km = $ab[0];
        $tshift_all = $ab[1];
        $pol_syn = $aa[37];
        $pol_obs = $aa[38];

        # test if weight or polarity exists. if neither then print nothing and dont skip space
        if (($aa[37]!=0  && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){

            $x = 0;
            @sensor_label = split('\.', $aa[0]);
            $inet = $sensor_label[1];
            $ista = $sensor_label[2];
            $iloc = $sensor_label[3];
            $icha = $sensor_label[4];
            # station label
            printf PLT "%f %f 10 0 0 1 $inet.$ista.$iloc.$icha\n", $x-0.8*$spis, $y;
            #printf PLT "%f %f 10 0 0 1 $sensor_label[1].$sensor_label[2].$sensor_label[3]\n", $x-0.8*$spis, $y;
            # station distance and overall time-shift
            if ($tshift_all==0.){
                printf PLT "%f %f 10 0 0 1 %d km\n", $x-0.8*$spis, $y-0.2, $dist_km;
                printf PLT "%f %f 10 0 0 1 %d\260 \n", $x-0.8*$spis, $y-0.4, $az[$i];  # azimuth (see az above)
            }
            else { 
                printf PLT "%f %f 10 0 0 1 %d km\n", $x-0.8*$spis, $y-0.2, $dist_km;
                printf PLT "%f %f 10 0 0 1 %d\260 \n", $x-0.8*$spis, $y-0.4, $az[$i];  # azimuth (see az above)
                printf PLT "%f %f 10 0 0 1 %.1f s\n", $x-0.8*$spis, $y-0.6, $tshift_all;  # tshift = Green_P_arrival - Input_P_arrival_weight_file
            }
            # azimuth
            # printf PLT "%f %f 10 0 0 1 %d\260 \n", $x-0.8*$spis, $y-0.6, $az[$i];  # azimuth (see az above)
            # polarities
            # NOTE if polarity is 0 or does not exist, then nothing is written
            if ($pol_syn || $keepBad==1) {
                if ($ab[1]==0.) {
                    printf PLT "%f %f 10 0 0 1 $pol_syn ($pol_obs)\n", $x-0.8*$spis, $y-0.6;
                }
                else {
                    printf PLT "%f %f 10 0 0 1 $pol_syn ($pol_obs)\n", $x-0.8*$spis, $y-0.8;
                }
            }
            $i=$i+1;
            $y--;
        } # end tests for weight and polarity
    }
    close(PLT);

    # plot data labels body waves
    open(PLT, $plt2_wf_info_b);
    $y = $nn-2;
    foreach (@capout_splice) {
      @aa = split;
      if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
      $x = 0;
      for($j=0;$j<2;$j+=$inc) {
          if ($aa[7*$j+2]>0 || $keepBad) {
            # printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;  # time shift each wf
            # printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;  # correl value
              $fracmis=sprintf("%2.2f", $aa[7*$j+3]);
              $logamp=sprintf("%2.2f", $aa[7*$j+6]);
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+5]\n", $x+0, $y-0.2;  # time shift each wf
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+4]\n", $x+0, $y-0.4;  # correl value
              printf PLT "%f %f 10 0 0 1 $fracmis\n", $x+0, $y-0.6;     # fractional misfit
              printf PLT "%f %f 10 0 0 1 $logamp\n", $x+0, $y-0.8;      # log(max_amp_data/max_amp_syn)
          }
          $x = $x + $x0[$j];
      }
      $y--;
    }
    # plot labels PR and PV
    $x = 0.2*$spib;
    for($j=0;$j<2;$j+=$inc) {
      printf PLT "%f %f 12 0 0 1 $name[$j]\n",$x,$nn-1.5;
      $x = $x+$x0[$j];
    }
  }
    close(PLT);

    # plot data labels surface waves
    open(PLT, $plt2_wf_info_s);
    $y = $nn-2;
    foreach (@capout_splice) {
      @aa = split;
      if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
#      $x = $x0[1];
      $x = 0;
      for($j=2;$j<5;$j+=$inc) {
          if ($aa[7*$j+2]>0 || $keepBad) {
              #printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;  # time shift each wave
              #printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;  # correl value
              $fracmis=sprintf("%2.2f", $aa[7*$j+3]);
              $logamp=sprintf("%2.2f", $aa[7*$j+6]);
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+5]\n", $x, $y-0.2;  # time shift each wave
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+4]\n", $x, $y-0.4;  # correl value
              printf PLT "%f %f 10 0 0 1 $fracmis\n", $x, $y-0.6;  # fractional misfit
              printf PLT "%f %f 10 0 0 1 $logamp\n", $x, $y-0.8;  # log(max_amp_data/max_amp_syn)
          }
          $x = $x + $x0[$j];     # original
      }
      $y--;
    }
    # -------------------- end plot data for each trace

    # plot labels PR PV SV SR SH for wave types

    # plot labels SV SR SH
    $x = 0.2*$spis;
#    $x = 0.2*$spis+$x0[2];
    for($j=2;$j<5;$j+=$inc) {
      printf PLT "%f %f 12 0 0 1 $name[$j]\n",$x,$nn-1.5;
      $x = $x+$x0[$j];
    }
  }
    close(PLT);

    
    #-----------------------------------------------------------
    print "Plotting little beachball ...\n";
    #-----------------------------------------------------------

    # note: magnitude scale is "fixed" at 1e17 for psmeca -Sm and 1 for psmeca -Sa
    open(PLT, $plt3);
    if ($tensor[1] eq "tensor") {
        # moment tensor is converted from AkiRichads basis to GCMT basis, which is required for psmeca
        printf PLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
    } else {
        # focal mechanism is plotted from the M0, strike/dip/rake values
        printf PLT "0 0 0 @meca[5,6,7] 1\n";  # 0.5*$spis,$nn-1;
    }
#    $x = 2;
#    foreach (@others) {
#       split;
#       printf PLT "%f -0.2 0 @_[1,2,3] 0.5 0 0 $_[6]\n",$x; $x+=1.5;
#    }
    close(PLT);

    #-----------------------------------------------------------
    print "Plotting station data, labels azimuths, weights, ...\n"; # (see staz above)
    #-----------------------------------------------------------
    #open(PLT, $plt4b);
    #foreach (@staz) {
    #  printf PLT;
    #}

    open(PLT, $plt4b);
    foreach (@staz_useweights) {
      printf PLT;
    }

    # Does this do anything??
    # plot station azimuths beachballs (see tkuh above)
    open(PLT, $plt4a);
    foreach (@tkuh) {
      printf PLT;
    }

    # plot piercing points on beachballs (see tklh above)
    #open(PLT, $plt4);
    #foreach (@tklh) {
    #  printf PLT;
    #}
    #close(PLT);

    open(PLT, $plt4);
    foreach (@tklh_useweights) {
      printf PLT;
    }
    close(PLT);

    # plot main label, 4 rows. next to the beachball.
    #$x = 0.5*$spis; 
    #$x = 2; 
    #$y = 0;
    #$tgap=0.5;
    ## plot four header labels (event type, focal mecha, var red, filters)
    ## Event 19910914190000000 Model 19910914190000000_wes_001 FM  350 56.985645  -74 Mw 5.80 rms 2.673e-06     1 CLVD -4.08 ISO  -4.464618 VR 7.8 data2 2.783e-06
    ##   0         1             2       3                     4    5      6        7  8   9   10      11      12  13    14  15       16    17  18   19     20
    #open(PLT, $plt4_5);
    #printf PLT "$x $y 12 0 0 Event $evid Model $model Depth $depth\n"; $y-=$tgap;
    #printf PLT "$x $y 12 0 0 @meca[4] %d %d %d @meca[8,9] @~g@~ %3.0f @~d@~ %3.0f @meca[10,11] VR %3.1f pol_wt %0.2f\n", @meca[5], @meca[6], @meca[7], @meca[14],@meca[16],@meca[18],$pol_wt;$y-=$tgap;
    #printf PLT "$x $y 12 0 0 $filterBand $duration\n" ; $y-=$tgap;  # 20120719 - filter bands
    #printf PLT "$x $y 12 0 0 @ncomp[1]" ;
    #close(PLT);
    #-----------------------------------------------------------
    # 2021-04-02 calvizuri - update for pstext gmt 6.1:
    ## Reads (x,y[,fontinfo,angle,justify],text) from <table> [or stdin].
    #$legend1="Event $evid Model $model Depth $depth\n";
    #$legend2="@meca[4] @meca[5], @meca[6], @meca[7],  @meca[8,9] @~g@~ @meca[14] @~d@~ @meca[16] @meca[10,11] VR @meca[18] pol_wt $pol_wt\n";
    #$legend3="$filterBand $duration\n" ;# $y-=$tgap;  # 20120719 - filter bands
    #$legend4="@ncomp[1]" ;
    ##printf PLT "1.0 1.0 Event $evid Model $model Depth $depth\n"; #$y-=$tgap;
    ##printf PLT "1.0 0.8 @meca[4] %d %d %d @meca[8,9] @~g@~ %3.0f @~d@~ %3.0f @meca[10,11] VR %3.1f pol_wt %0.2f\n",
    ##                    @meca[5], @meca[6], @meca[7], @meca[14],@meca[16],@meca[18],$pol_wt; #$y-=$tgap;
    ##printf PLT "1.0 0.6 $filterBand $duration\n" ;# $y-=$tgap;  # 20120719 - filter bands
    ##printf PLT "1.0 0.4 @ncomp[1]" ;
    #open(PLT, $plt4_5);
    #printf PLT "10 1 12 0 0 0 $legend1";
    #printf PLT "10 1 12 0 0 0 $legend2";
    #printf PLT "10 1 12 0 0 0 $legend3";
    #printf PLT "10 1 12 0 0 0 $legend4";
    #close(PLT);
    #printf STDERR "$x $y 12 0 0 $filterBand $duration y $y tgap $tgap \n";  # 20120719 - filter bands
    #-----------------------------------------------------------
    ## TRY 3. works. Original code messy. Needed: -F+jl to justify, plus x y. M flag not needed. weird.
    print "Plotting header labels ...\n";
    open(PLT, $plt4_5);
    # > 0 -0.5 13p 3i l
    printf PLT "0.0 -0.0 Event $evid Model $model Depth $depth\n";
    printf PLT "0.0 -0.5 @meca[4] %d %d %d @meca[8,9] @~g@~ %3.0f @~d@~ %3.0f @meca[10,11] VR %3.1f pol_wt %0.2f\n", @meca[5], @meca[6], @meca[7], @meca[14],@meca[16],@meca[18],$pol_wt;;
    printf PLT "0.0 -1.0 $filterBand $duration\n" ; # 20120719 - filter bands
    printf PLT "0.0 -1.5 @ncomp[1]" ;
    close(PLT);

    #-----------------------------------------------------------

  print "Done.\n";
  print "============================\n";
  }  # while (@capout) {

  #---------------------------------
  # FIGURE 2: big moment tensor with station names at lower-hemisphere piercing points
  print STDERR "Plotting big beachball ...\n";
  #---------------------------------
  $pwidth_in = 8.5;  # width of paper
  $pheight_in = 11;  # height of paper
  #system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");
  system("gmt gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");

  # restore
  @capout=@capout0;

  while (@capout) {
    @capout_splice = splice(@capout,0,$nn-2);

    # plot beachball (see notes above)
    open(XPLT, $xplt3);
    if ($tensor[1] eq "tensor") {
     printf XPLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
    } else {
     printf XPLT "0 0 0 @meca[5,6,7] 1\n"; #0.5*$spis,$nn-1;
    }
    close(XPLT);

    # plot piercing points on beachballs (see tklh above)
    $i=0; $j=0; $k=0;
    open(XPLT, $xplt4);
    open(XPLTC, $xplt4c);
    open(XPLTD, $xplt4d);
    open(XPLTE, $xplt4e);
    open(XPLTF, $xplt4f);
    foreach (@tklh_useweights) {
	if ($ifmp[$i] * $ifmpt[$i] < 0) {     # mismatcing polarities
	    if ($ifmp[$i]>0){printf XPLTE;}   # input is UP (+1); theoretical is DOWN (-1)
	    else {printf XPLTF;}}             # input is DOWN (-1); theoretical is UP (+1)
	elsif ($ifmp[$i]>0){printf XPLTC;}    # both input and theoretical are UP (+1)
	elsif ($ifmp[$i]<0){printf XPLTD;}    # both input and theoretical are DOWN (-1)
        else {printf XPLT;}                   # no input polarity pick in the weight file
	$i=$i+1;
    }
    close(XPLT);
    close(XPLTC);
    close(XPLTD);
    close(XPLTE);
    close(XPLTF);

    print "Plotting azimuths and station name ...\n";
if ($only_pol == 0) {
    # plot station azimuths beachballs (see staz above)
    open(XPLT, $xplt4b);
    foreach (@staz_useweights) {
      printf XPLT;
    }
    close(XPLT);

    # plot station azimuths beachballs (see tkuh above)
    open(XPLT, $xplt4a);
    foreach (@tkuh_useweights) {
      printf XPLT;
    }
    close(XPLT);
#------------

    open(XPLT, $xplt5a);
    foreach (@staz_useweights) {
      @aa = split;
      @aa_split = split('\.', $aa[2]);
      printf XPLT "%s %s 8 0 0 CB %s.%s.%s\n", 
          $aa[0], $aa[1], $aa_split[1], $aa_split[2], $aa_split[3];
    }
    close(XPLT);

#     open(XPLT, $xplt5b);
#     foreach (@tkuh) {
#       @aa = split;
#       printf XPLT "%s %s 8 0 0 CB (%s)\n",$aa[0],$aa[1],$aa[2]; 
#     }
#     close(XPLT);

    open(XPLT, $xplt5c);
    foreach (@tklh_useweights) {
      @aa = split;
      @aa_split = split('\.', $aa[2]);
      printf XPLT "%s %s 8 0 0 CB %s.%s.%s\n", 
          $aa[0], $aa[1], $aa_split[1], $aa_split[2], $aa_split[3];
    }
    close(XPLT);

    print "Plotting title ...\n";
    $x = -1; 
    $y = 0;
    open(XPLT, $xplt6);
    printf XPLT "0 0 16 0 0 0 @meca[0..3]\n";
    # Event 19910914190000000 Model 19910914190000000_wes_001 FM  350 56.985645  -74 Mw 5.80 rms 2.673e-06     1 CLVD -4.08 ISO  -4.464618 VR 7.8 data2 2.783e-06
    #   0         1             2       3                     4    5      6        7  8   9   10      11      12  13    14  15       16    17  18   19     20
    printf XPLT "0 -0.05 16 0 0 0 @meca[4] %d %d %d @meca[8,9] @~g@~ %3.0f @~d@~ %3.0f @meca[10,11] VR %3.1f pol_wt %0.2f\n",@meca[5], @meca[6], @meca[7], @meca[14],@meca[16],@meca[18], $pol_wt;
    close(XPLT);

}
  print STDERR "Done.\n";
  print "============================\n";
}
#---------------------------------
  print "Summary results: @meca\n\n";

}
1;
