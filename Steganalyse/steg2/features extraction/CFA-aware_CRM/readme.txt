A few notes to CFA-aware feature sets as proposed in [2] as a modification of SRMQ1.
------------------------------------------------------------------------------------

Input image must have "blue" pixels (see Bayer CFA) in the first row, first column. 
NII/INI ..... SRMQ1CFA.m
RB/GG  ...... SRMQ1C2cfa.m
R/B/GG ...... SRMQ1C2cfaRGGB.m


SRMQ1CFAred.m is a modified version of SRMQ1CFA.m allowing to input any color image with explicitely identified "red" pixel in the fourth input of SRMQ1CFAred.
f = SRMQ1CFAred(IMAGE,Tc,part,4) is the same as f = SRMQ1CFA(IMAGE,Tc,part)

The input parameter part must be set to 'color', Tc = 2.

Example:
X = double(imread('image_name.ppm'));
Features = SRMQ1CFA(X,2,'color');	% feature structure
F = struct2array(Features);    		% feature vector (length(F)==5514)


"Magnificient 8" are 8 features (out of 200) in Features.s1_minmax41c_q1 with indices [94,95,99,100,119,120,124,125], i.e. cooccurrence bins 
   1     2     1     2     1     2     1     2 
   1     1     2     2     1     1     2     2 
   1     1     1     1     2     2     2     2 
"Magnificient 8" can also be extracted from F as F([94,95,99,100,119,120,124,125]+702).


Reference:
[2] CFA-aware Features for Steganalysis of Color Images, M. Goljan and 
J. Fridrich, Proceedings SPIE, Electronic Imaging, Media Forensics
and Security, 2015.
