function fout = rationaliseModFreq(fin)
% rationaliseModFreq
%
% Most of the dataset are 50, 150Hz (100Hz) steps. 
% Rationalised frequency shifts them all to this for 
% ease of population analysis. 
% This mainly affects the minority of data which are 
% in even multiples of 100Hz. These are rounded down to the 
% nearest x50Hz. 

fout = 100*ceil(fin/100)-50;