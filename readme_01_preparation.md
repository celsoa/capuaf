# FMTU: full moment tensor and uncertainty estimation using seismic waveforms

## Earthquake waveform data 
The waveform data should be prepared using the python utility ``pysep``, which will fetch the waveforms from repositories, perform various tests and prepare the data for analysis, including correcting for instrument response, rotating to source-receiver coordinate system, and preparing auxiliary files used in the moment tensor estimation.

All waveforms will be saved in an event directory named after the origin time (e.g. ``20100516063454464``).
The waveforms should be in units of velocity in cm/s or displacement in cm, and in SAC format.
The SAC headers should include a reference time set as the origin time, epicentral distance and azimuth. 

## Weight files
There should be a file named weight.dat in the same directory, in the following format:
```
station_name    dist    w1 w2 w3 w4 w5  tp tp_w ts ts_w ti_r [ti_s]
```
Here dist specifies the source-receiver distance in km, and is matched with the corresponding Greens functions (dist.grn.?).
If first motion polarities are available, they need to be entered directly after station_name, followed by a `/` and value of the polarity (up=1, down=-1).
The weights w1 to w5 set the relative importance of each of the 5 waveform windows during moment tensor estimation, see ``Time window determination`` below for details.
The window times `txx` are:
```
tp:   time of the first P arrival, if it's set to a positive value
tp_w: length of the body wave window
ts:   arrival time for surface waves
ts_w: length of the surface wave window
ti_r, ti_s: are the initial time shifts for the surface waves (rayleigh, love)
``` 
Positive shifts mean the observed is delayed with respect to the model.
If w2 is set to -1, it indicates that the station is at teleseismic distances and only the P (PnlZ) and SH (T) are used.
In this case, ts is the S arrival time when it is positive.

## Greens functions library
The Greens functions are computed using the FK routines and named (automatically) in the format ``dist.grn.[a-c,0-8]``, where dist is the source-receiver distance in km and [a-c, 0-8] are the 9 Greens functions components.
All Greens functions from one source depth are placed in a single directory named in the format ``model_depth``.
The Greens functions are created in SAC format with two time marks set: t1 for the first P arrival and t2 for the first S arrival.
If first-motion data are to be used in the inversion, the Greens functions need to have user1 and user2 set as the P and S take-off angles (in degrees from down).

## Time window determination
The moment tensor estimation routines work by cutting each station waveform into windows and aligning them individually with the synthetic waveforms.
Depending on data availability, the waveforms for each station will be cut-up in up to 5 waveform windows: PnlZ, PnlR, Z, R, T.
These are 2 Pnl windows (vertical and radial components), and 3 surface waveform windows (Surf-vertical, -radial, -transverse; or Rayleigh vertical and horizontal, and one transverse).
The windows are determined in one of the following ways:
1) If the SAC headers have t1 and t2 time marks, the code will use them for determining the Pnl window. The same is true for the surface wave window (using t3 and t4 headers).
2) If positive apparent velocities are given in the code (see flag in cap.pl), they will used to calculate the time windows as follows:
```
t1 = dist/vp        - 0.3*tmax_body
t2 = ts             + 0.2*tmax_body
t3 = dist/vLove     - 0.3*tmax_surf
t4 = dist/vRayleigh + 0.7*tmax_surf
```
3) Using the SAC headers tp, ts in the Greens functions
```
t1 = tp - 0.2*tmax_body
t2 = t1 +     tmax_body
t3 = ts - 0.3*tmax_surf
t4 = t3 +     tmax_surf
```
Here tmax_body, tmax_surf are the maximum lengths for the Pnl and surface waves windows
(see the -T options below).
