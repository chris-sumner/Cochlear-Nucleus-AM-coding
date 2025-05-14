# Cochlear-Nucleus-AM-coding - PLOS Biology version
The data and code for envelope coding in the cochlear nucleus. This requires MATLAB. It also contains some C (Mex) code which has been compiled for 64-bit Windows (so some work required to run elsewhere). 

This is ALL the data and code to go from spike times to published figures. Note that some of the analyses take a long time (SAC peak picking took weeks to run) to run from scratch, but processed data files will allow reproduction of most of the figures quickly. 

## Folder structure
`Root` - contains the .m files from which eveyrthing is run, from analysis to figure production.

`rawdata` - the raw data. Spike times for AM data, information about the data, and the example pure tone responses used. 

`functions` - a set of supporting functions. May contain some that are not used. 

`datafiles` - the data in processed form - as produced by the code below. Not necessary as you can remake. Note that some of these file were too large to fit on Github in original form and here we provide minimal versions of these. One intermediate file (which was ~1.3GB) and is not provided, but it not required to reproduce the figures. 

`figure_datavalues` - Excel spreadhsheets of the raw data values underlying any figure panels which show summary values (means etc.), as requested by PLOS Biology.  

## Main code files in root folder
There are two files which peform the core data analysis:

`SpikeStatsCalc.m`  - this loads up the data and computes various statistics such as the Vector Strength, mode-locking analysis (Z-ISI) and picking of SAC peaks. 

`WolgRonacClassification.m` - this performs the spike train classifier.

Both write .mat files to datafiles\ for use in subsequent analysis/display.

The remaining files pull together, do additional analyses, and make the figures:

`ExampleUnitFigures.m` - this makes the two example figures, and the methods figure. 

`mainAnalyses.m` - integrates the above analyses files, performs the additional analyses and makes the remaining figures. Start here. 

`populationCoding.m` - this runs the classifier for small clusters of neurons. Run from within mainAnalyses.

`groupSpikeStatsByMTF.m` - wrangles the data structures from SpikeStatsCalc in to a form where each data item is a single modulation transfer function. 

## Raw data
The raw data folder contains the original data as follows.

The data come from two previous studies (Rhode and Greenberg 1994, Rhode 1994) and were collected in 1988 and 1991 respectively. They contain all the data from the original papers, including some data we have not analysed in our study (e.g. responses to AM with a broadband noise carrier, >100% modulation, and neurons with low characteristic frequencies).

Spike times in response to AM tones are stored in four separate .mat files (`justSpikeTimesA`...`D`), in order to meet the Github file size limit (<25MB). The function loadSpikeTimes loads these and joins them together correctly for all data analysis. Each entry (in the cell array) is a response of a single neuron to a single stimulus condition. Each repeat is a different sub-cell. There are ~57000 entries. 

The information about the dataset is stored in `unitList.mat`. It contains the following variables:

`wholeList171212` - a large matrix, with a row for each entry in the above array, which describes each dataset.

`EXPLOGLIST` - is a cell array with a text descriptor for each column in this. The most important of these are:

- **Expt** - subject ID. 
- **Unit** - the number of the neuron within the Expt.
- **typeNum** - a numeric code to the pure tone classification for the neuron. 
- **cf** - the characteristic frequency of the neuron.
- **cf_thr** - the threshold sond lelve of the neuron at CF.
- **mod_level** - the sound level (dB SPL) of the AM stimulus.
- **mod_freq** - the modulation frequency of the AM stimulus (in Hz).
- **amToneDur** - the duration (ms) of the AM stimulus.
- **carrFreq** - the carrier frequency of the tone (Hz).
- **depthMod** - modulation depth of the AM stimulus. 
  
There are other statistics but these have not been checked recently so are not described here. 

`typeNames` - the pure tone unit types as abbreviated text, which TypeNum correspond to. 
  
Rhode, W.S. and S. Greenberg. J Neurophysiol, 1994. 71(5): p. 1797-825.

Rhode, W.S. Hear Res, 1994. 77(1-2): p. 43-68.




