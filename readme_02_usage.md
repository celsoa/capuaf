# FMTU: full moment tensor and uncertainty estimation using seismic waveforms

## Examples

The following commands find the best focal mechanism and moment magnitude for the 2008/4/18 Southern Illinois earthquake 20080418093700 using the central US crustal velocity model cus with the earthquake at a depth of 15 km.
The examples assume that the Greens functions have already been computed and saved in $green/cus/cus_15/.

### EXAMPLE: RANDOM SEARCH
* Double Couple
``` 
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I100000 20080418093700
```

* Full Moment Tensor
```
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1       -Y1 -I100000 20080418093700
```

### EXAMPLE: GRID SEARCH
* Double Couple
```
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I1/1/37/10/19 20080418093700
```
* Full Moment Tensor
```
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1       -Y1 -I10/10/37/10/19 20080418093700
```

### EXAMPLE: DEPTH SEARCH
To find the best focal depth, repeat the inversion for different focal depths either by using a for-loop or the '-A' flag
* for loop with depths 5 to 30 km at 5-km intervals
```
> for h in 5 10 15 20 25 30; do cap.pl -H0.2 -P1 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_$h/5.0 -E0 -K0 -Y2 20080418093700; done
```
* same as above but using the A flag in the format -Astart/end/incr
```
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I1/1/37/10/19 20080418093700 -A5/30/5
```

### EXAMPLE: PLOT THE DEPTH SEARCH RESULTS
```
> depth_test 20080418093700 cus
> gv dep_20080418093700.ps
```

## RESULTS
The examples above will produce postscript figures with waveform fits and
beachball, and an output file (extension .out). The waveform fits are saved in
the event directory, in file cus_15.ps (double-couple) or cus_15_fmt.ps (full
moment tensor).

### OUTPUT FILE AND HEADER INFO
The results are saved in cus_15.out, which has a header like the following:
```
Event 20080418093700 Model cus_015 FM  297 86.815262    0 Mw 5.20 rms 4.547e-05   112 ERR   0   0   0 CLVD -1.89 -nan ISO  10.212961 0.00 VR 80.3 data2 1.026e-04
# Hypocenter_sac_header elat 3.845000e+01 elon -8.789000e+01 edep 1.160000e+01
# tensor = 7.943e+23  0.9511 -0.5865 -0.0273 -0.6243  0.0474  0.1075
# norm L1    # Pwin 35 Swin 70    # N 8 Np 16 Ns 24
```

### HEADER INFO
```
H1: event ID, model, depth. FM: focal mechanism with strike 297, dip 87, rake 0.
H1: The values CLVD -1.89, ISO  10.212961 correspond to lune coordinates (gamma, delta) in degrees (if resolving a full moment tensor; for DC they are zero).
H1: rms is misfit, Data2 is data norm, VR is Variance reduction.
H2: Hypocenter location, obtained from SAC headers.
H3: Moment tensor, format: Mxx Mxy Mxz Myy Myz Mzz (where x=North, y=East, z=Down).
H4: Norm (L1 or L2); Window length of P (Pwin) and S (surf); Number of stations (N); Number of P windows (Np) and Number of Surf windows (Ns).
The rest of the output includes individual station data, eg rms, cross-correlation coef., time shift of individual waveforms, etc. 
For a detailed description see Alvizuri et al. (2018).
```
