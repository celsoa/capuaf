# _FMTU_: a fast seismic source estimation package 
* fetch worldwide earthquake waveforms
* estimate seismic full moment tensors
* estimate moment tensor uncertainties
* analyze source types

## Dependencies 

### Waveform fetching and preparation
Earthquake waveform data can be obtained with the ``pysep`` utility, which is a python package for fetching waveform data through FDSN and EIDA webservices, and preparing it for moment tensor estimation.
Thepackage can be obtained from the github repository:
```
git clone XXX
```

### Signal processing libraries 
* ``libsac.c``: this is a library of routines for filtering waveform data.
* The library is distributed with the Seismic Analysis Code (SAC), which can be requested to IRIS at: https://ds.iris.edu/ds/nodes/dmc/forms/sac/

### Greens functions
* CAP requires a library of precomputed greens functions which are computed with the ``fk`` routines
* The latest version ``fk3.3`` can be requested at http://www.eas.slu.edu/People/LZhu/home.html 

## Installation
### Compilation options
```
make cap d=WB               # Compile cap to write a binary file with search data
make cap f=omp              # Compile cap to run in parallel mode (openMP)
make cap d=WB f=omp         # Compile cap to run in parallel and write a binary datafile 
make clean cap d=WB f=omp   # Clean all and compile with all flags
```

### Paths
* executables
* waveform data

## Usage
* A typical command to estimate an earthquake  mechanism looks like this
```
cap.pl EVID [flag1/varA] [flag2/varB] ... [flagN/varX]
```
* Many flags and values can be passed to FMTU. Which flags and their values depend on the earthquake type, the earthquake signal, etc.
* The following example command will estimate the focal mechanism for the 2008-04-18 earthquake in central USA
```
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -Y1 -I100000 20080418093700
```
* For more details and to get started:
    - type ``cap.pl`` in the command line and see the flags there
    - see [readme_preparation.md](#readme_preparation.md) and [readme_usage.md](#readme_usage.md)

## References
If you use this code in your research, please cite the first or all of the following
``` 
Estimation of full moment tensors, including uncertainties, for nuclear explosions, 
volcanic events, and earthquakes. Journal of Geophysical Research: Solid Earth, 
123(6), 5099-5119. Alvizuri, C., Silwal, V., Krischer, L., & Tape,C. (2018).

Alvizuri, C., & Tape, C. (2016). Full moment tensors for small events (Mw<3)
at Uturuncu volcano, Bolivia. Geophysical Journal International, 206(3),
1761-1783.

Silwal, V., & Tape, C. (2016). Seismic moment tensors and estimated
uncertainties in southern Alaska. Journal of Geophysical Research: Solid Earth,
121(4), 2772-2797.

These may also be appropriate
Zhu and Ben-Zion, 2013, GJI
Zhu and Helmberger, 1996, BSSA 86, 1645-1641
Zhao and Helmberger, 1994, BSSA 84(1), 91-104
```

