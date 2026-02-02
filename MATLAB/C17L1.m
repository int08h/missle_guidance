clear
global a gm count
global ArrayICOUNT ArrayBGAM ArrayVRX ArrayVRY ArrayTF 
format long e
XLONGMDEG=45.;
XLONGTDEG=90.;
ALTNMT=0.;
ALTNMM=0.;
TF=1000;
DEGRAD=360./(2.*pi);
a=2.0926e7;
gm=1.4077e16;
ALTT=ALTNMT*6076.;
ALTM=ALTNMM*6076.;
XLONGM=XLONGMDEG/DEGRAD;
XLONGT=XLONGTDEG/DEGRAD;
XM=(a+ALTM)*cos(XLONGM);
YM=(a+ALTM)*sin(XLONGM);
XT=(a+ALTT)*cos(XLONGT);
YT=(a+ALTT)*sin(XLONGT);
[VRXM,VRYM]=olambert(XM,YM,TF,XT,YT,XLONGM,XLONGT)
 output=[ArrayVRX', ArrayVRY', ArrayTF'];
 disp('The final iteration')
 count
VRXM
VRYM
% data to file
 save datfil.txt output -ascii;
