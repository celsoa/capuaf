# this subroutine plots waveform fits produced by source inversion srct

sub plot {

#  local($mdl, $t1, $t2, $am, $num_com, $sec_per_inch) = @_; # original
  local($mdl, $t1, $t2, $am, $num_com, $sec_per_inch, $filterBand, $fmt_flag) = @_;
  local($nn,$tt,$plt1,$plt2,$plt3,$plt4,$i,$nam,$com1,$com2,$j,$x,$y,@aa,$rslt,@name,@aztk);
  local $keepBad = 0;
  
  @trace = ("1/255/255/255","3/0/0/0");       # plot data trace
  @name = ("P V","P R","Surf. V"," Surf R","SH");

  $filterBand = "Filter periods (seconds): $filterBand";    # 20120719 - report filter bands
  
#--------------------------
# plotting GMT file with synthetics and data
  
  # set GMT defaults
  # for all options (such as PAPER_MEDIA): http://gmt.soest.hawaii.edu/gmt/html/man/gmtdefaults.html
  # to check defaults, e.g.: gmtdefaults -L | grep MEASURE_UNIT
  @dum = split('_', $mdl);  # split mdl string
  $outfile = sprintf("%s_%03d.out",@dum[0],int(@dum[1]));

  # read in the output file results
#  open(FFF,"$mdl.out"); # original
  open(FFF,$outfile);    # 20130102 calvizuri - new file name
  @rslt = <FFF>;
  close(FFF);
  @meca = split('\s+',shift(@rslt));
  @variance = split('\s+',shift(@rslt));
  @tensor = split('\s+',$rslt[0]);
  @others = grep(/^#/,@rslt);
  @rslt=grep(!/^#/,@rslt);
  $nrow = @rslt;

  $pwidth_in = 8.5;  # width of paper
  $pheight_in = $nrow + 2;  # height of paper
  print "\n$nrow rows to plot";
  print "\npaper is $pwidth_in inches wide and $pheight_in inches tall";
  system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");

  # positions of seismograms on the page

  # height of each seismogram
  $nn = int($pheight_in);
  $height = $pheight_in - 0.5;
  #($nn,$height) = (12,10.5);   # 10 rows of traces per 10.5 in.
  #print "\n$nn rows of traces per $height in";  
  print "\nseconds per inch = $sec_per_inch";
  $sepa = 0.1*$sec_per_inch;
  ($tt, $inc) = (2*$t1 + 3*$t2 + 4*$sepa, 1);
  ($tt, $inc) = ($t1 + $t2 + $sepa, 4) if $num_com == 2;
  $width = 0.1*int(10*$tt/$sec_per_inch+0.5);
  @x0 = ($t1+$sepa, $t1+$sepa, $t2+$sepa, $t2+$sepa, $t2);

# pssac2 amplitude scaling:
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

  # KEY: set amplitude scaling for seismograms
  if ($am>0.) {$stam = "$am/-1";} else {$stam=-$am;} # original line (with pssac, not pssac2)
  $stam = "$am/0";                                   # overwrite for absolute amplitudes (to match default plotting)
  print "\namplitude scaling am = $am";
  print "\npssac2 amplitude scaling stam = $stam\n";
#  $outps = "$mdl.ps";   # original
  $outps = sprintf("%s_%03d.ps",@dum[0],int(@dum[1])); # reformatted filename
  $outps = sprintf("%s_%03d_fmt.ps",@dum[0],int(@dum[1])) if $fmt_flag eq "true";

  # (1) plot cut seismograms with scaled amplitudes (first command: no -O appears)
  $plt1 = "| pssac2 -JX${width}i/${height}i -R0/$tt/0/$nn -Y0.2i -Ent-2 -M$stam -K -P >> $outps";

  # (2) plot text labels
  $plt2 = "| pstext -JX -R -O -K -N >> $outps";

  # (3) plot beachballs (solution, followed by possible local minima)
  $ballcolor = "150";
  $dY = ${pheight_in} - 1.8;
  $plt3 = "| psmeca -JX5i/1i -R-1/9/-1/1 -Sa5i -G$ballcolor -Y${dY}i -X-0.7i -O -K >> $outps";
  $plt3 = "| psmeca -JX5i/1i -R-1/9/-1/1 -Sm8i -G$ballcolor -Y${dY}i -X-0.7i -O -K >> $outps" if $tensor[1] eq "tensor";

  # (4) plot markers on beachball
  # note: -JPa is a basemap for polar coordinates, clockwise from north

  # azimuths
  $plt4b = "| psxy -JPa1i -R0/360/0/1 -Sc0.02i -N -W0.5p,0/0/0 -G255 -O -K >> $outps";

  # supplemental: upper hemisphere piercing points on beachballs (o)
  $plt4a = "| psxy -JPa1i -R0/360/0/1 -Sc0.08i -N -W0.5p,255/0/0 -O -K >> $outps";

  # default: lower hemisphere piercing points on beachballs (x) (last command: no -K appears)
  $plt4 = "| psxy -JPa1i -R0/360/0/1 -Sx0.10i -N -W0.5p,255/0/0 -G255 -O >> $outps";

  #$plt1=$plt2=$plt3="|cat";	# output GMT commands to command window for testing

#--------------------------

#  $outps2 = "${mdl}_beach.ps"; # original
  $outps2 = sprintf("%s_%03d_beach.ps",@dum[0],int(@dum[1]));   # 20130102 calvizuri - revised filename
  $outps2 = sprintf("%s_%03d_beach_fmt.ps",@dum[0],int(@dum[1])) if $fmt_flag eq "true";

  $fac = 6.5;
  $fac2 = 8.2*$fac;   # original: 5*$fac
  $JP = "-JPa${fac}i";

  # plot beachball
# $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sm${fac2}i -X1i -Y2i -K -P >> $outps2";
  $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sa${fac}i -X1i -Y2i -K -P >> $outps2";
  $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sm${fac2}i -X1i -Y2i -K -P >> $outps2" if $tensor[1] eq "tensor";

  # plot markers on beachball
  # note: -JPa is a basemap for polar coordinates, clockwise from north

  # azimuths
  $xplt4b = "| psxy $JP -R0/360/0/1 -Sc0.02i -N -W0.5p,0/0/0 -G255 -O -K >> $outps2";

  # supplemental: upper hemisphere piercing points on beachballs (o)
  $xplt4a = "| psxy $JP -R0/360/0/1 -Sc0.08i -N -W0.5p,255/0/0 -O -K >> $outps2";

  # default: lower hemisphere piercing points on beachballs (x)
  $xplt4 = "| psxy $JP -R0/360/0/1 -Sx0.10i -N -W0.5p,255/0/0 -G255 -O -K >> $outps2";

  # plot text labels
  $xplt5a = "| pstext $JP -R0/360/0/1 -N -O -K >> $outps2";
  $xplt5b = "| pstext $JP -R0/360/0/1 -N -O -K >> $outps2";
  $xplt5c = "| pstext $JP -R0/360/0/1 -N -O -K >> $outps2";

  # title (LAST COMMAND: no -K appears)
  $xplt6 = "| pstext -JX -R -N -O -Xa0 -Ya7.5 >> $outps2";

#--------------------------
# FIGURE 1: waveform fits with moment tensor

  # uncomment output for debugging purposes
  print "\n-------------------";
  print "\nmeca:\n@meca";
  print "\nvariance:\n@variance";
  print "\ntensor:\n@tensor";
  print "\nothers:\n@others";
  print "\nrslt:\n@rslt";
  print "\n-------------------\n";

  # compute piercing points for beachballs
  $i = 0; $j = 0;
  $pi = 3.14159265358979323846;
  @tklh=(); @tkuh=(); @staz=(); @az=();
  foreach (@rslt) {
    @aa = split;
    $stnm = $aa[0];                              # station name
    #next if $aa[2] == 0;                        # skip if no body waves
    $x = `saclst az user1 f ${mdl}_$aa[0].0`;    # get the azimuth and P take-off angle
    @aa = split(' ', $x);
    #print "\n--> saclst az user1 f ${mdl}_$aa[0].0";

    # compute polar coordinates azimuth and radius
    $az[$i] = $aa[1];
    $staz[$i] = sprintf("%s %f %s\n",$aa[1],1.1,$stnm);      # station azimuth
    if ($aa[2]>90.) {                                         # upper hemisphere
       $rad = sqrt(2.)*cos($aa[2]*$pi/360);
       $tkuh[$j] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);
       $j++;

       # project piercing point to lower hemisphere
       $aa[1] += 180;
       $aa[2]=180-$aa[2];
    }
    $rad = sqrt(2.)*sin($aa[2]*$pi/360);
    $tklh[$i] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);        # lower hemisphere
    $i++;
  }

  # remove the file if it exists
  unlink($outps) if -e $outps;
  unlink($outps2) if -e $outps2;

  # save a copy for the second file
  @rslt0=@rslt;

  while (@rslt) {

    # plot waveforms
    open(PLT, $plt1);
    $i = 0;
    @aaaa = splice(@rslt,0,$nn-2);
    foreach (@aaaa) {
      @aa = split;
      $nam = "${mdl}_$aa[0].";
      $x=0;
      for($j=0;$j<5;$j+=$inc) {
        $com1=8-2*$j; $com2=$com1+1;
	if ($aa[4*$j+2]>0) {
	   printf PLT "%s %f %f 5/0/0/0\n",$nam.$com1,$x,$nn-$i-2;
           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
	} elsif ($keepBad) {
	   printf PLT "%s %f %f 2/0/255/0\n",$nam.$com1,$x,$nn-$i-2;
           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
	}
        $x = $x + $x0[$j];
      }
      $i++;
    }
    close(PLT);
    
    # text labels
    open(PLT, $plt2);
    $y = $nn-2;
    $i=0;
    foreach (@aaaa) {
      @aa = split;
      #print "\n @aa";
      $x = 0;
      printf PLT "%f %f 10 0 0 1 $aa[0]\n",$x-0.8*$sec_per_inch,$y;            # station label
      printf PLT "%f %f 10 0 0 1 $aa[1]\n",$x-0.7*$sec_per_inch,$y-0.2;        # distance_km/overal time shift
      printf PLT "%f %f 10 0 0 1 %.1f\n",$x-0.7*$sec_per_inch,$y-0.4,$az[$i];  # azimuth (see az above)
      $i=$i+1;
      for($j=0;$j<5;$j+=$inc) {
	if ($aa[4*$j+2]>0 || $keepBad) {
          printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;
          printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;
	}
        $x = $x + $x0[$j];
      }
      $y--;
    }
    $x = 0.5*$sec_per_inch; 
    $y = $nn-0.2;
    printf PLT "$x $y 12 0 0 0 @meca[0,1,2] and Depth $meca[3]\n"; $y-=0.3;
    printf PLT "$x $y 12 0 0 0 @meca[4..22]\n";$y-=0.3;
    printf PLT "$x $y 12 0 0 0 @variance[1..3]\n" if $variance[1] eq "Variance";
#    printf PLT "%f %f 10 0 0 0 @meca\n",0.5*$sec_per_inch,$nn-0.4;  # full title    # original
#    printf PLT "%f %f 10 0 0 0 @meca\n",0.5*$sec_per_inch,$nn-0.2;  # full title
    printf PLT "%f %f 12 0 0 0 $filterBand.\n",0.5*$sec_per_inch,$nn-1.1;  # 20120719 - filter bands
    $x = 0.2*$sec_per_inch;
    for($j=0;$j<5;$j+=$inc) {
      printf PLT "%f %f 12 0 0 1 $name[$j]\n",$x,$nn-1.5;
      $x = $x+$x0[$j];
    }
    close(PLT);
    
    # plot beachball
    # note 1: magnitude scale is "fixed" at 1e17 for psmeca -Sm and 1 for psmeca -Sa
    # note 2: moment tensor is converted from AkiRichads basis to GCMT basis, which is required for psmeca
    open(PLT, $plt3);
    if ($tensor[1] eq "tensor") {
        printf PLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
    } else {
        printf PLT "0 0 0 @meca[5,6,7] 1\n";  # 0.5*$sec_per_inch,$nn-1;
    }
#    $x = 2;
#    foreach (@others) {
#       split;
#       printf PLT "%f -0.2 0 @_[1,2,3] 0.5 0 0 $_[6]\n",$x; $x+=1.5;
#    }
    close(PLT);

    # plot station azimuths beachballs (see staz above)
    open(PLT, $plt4b);
    foreach (@staz) {
      printf PLT;
    }

    # plot station azimuths beachballs (see tkuh above)
    open(PLT, $plt4a);
    foreach (@tkuh) {
      printf PLT;
    }

    # plot piercing points on beachballs (see tklh above)
    open(PLT, $plt4);
    foreach (@tklh) {
      printf PLT;
    }
    close(PLT);
    
  }  # while (@rslt) {

#---------------------------------
# FIGURE 2: big moment tensor with station names at lower-hemisphere piercing points

  $pwidth_in = 8.5;  # width of paper
  $pheight_in = 11;  # height of paper
  system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");

  # restore
  @rslt=@rslt0;

  while (@rslt) {
    @aaaa = splice(@rslt,0,$nn-2);
    
    # plot beachball (see notes above)
    open(XPLT, $xplt3);
    if ($tensor[1] eq "tensor") {
     printf XPLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
    } else {
     printf XPLT "0 0 0 @meca[5,6,7] 1\n"; #0.5*$sec_per_inch,$nn-1;
    }
    close(XPLT);

    # plot station azimuths beachballs (see staz above)
    open(XPLT, $xplt4b);
    foreach (@staz) {
      printf XPLT;
    }

    # plot station azimuths beachballs (see tkuh above)
    open(XPLT, $xplt4a);
    foreach (@tkuh) {
      printf XPLT;
    }

    # plot piercing points on beachballs (see tklh above)
    open(XPLT, $xplt4);
    foreach (@tklh) {
      printf XPLT;
    }
    close(XPLT);

#------------

    open(XPLT, $xplt5a);
    foreach (@staz) {
      @aa = split;
      printf XPLT "%s %s 8 0 0 CB %s\n",$aa[0],$aa[1],$aa[2]; 
    }
    close(XPLT);

#     open(XPLT, $xplt5b);
#     foreach (@tkuh) {
#       @aa = split;
#       printf XPLT "%s %s 8 0 0 CB (%s)\n",$aa[0],$aa[1],$aa[2]; 
#     }
#     close(XPLT);

    open(XPLT, $xplt5c);
    foreach (@tklh) {
      @aa = split;
      printf XPLT "%s %s 8 0 0 CB %s\n",$aa[0],$aa[1],$aa[2]; 
    }
    close(XPLT);

    # TITLE
    open(XPLT, $xplt6);
    printf XPLT "0 0 18 0 0 0 @meca[0..9]\n",0,0;
    printf XPLT "0 -0.05 18 0 0 0 @meca[10..22]\n",0,0;
    close(XPLT);

  }  # while (@rslt) {

#---------------------------------

}
1;
