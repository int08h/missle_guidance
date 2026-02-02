clear
count=0;
for I=2:160
W=10^(.025*I-1);
C1=1.;
C2=.0363;
TOP=9000*sqrt(C1^2+(C2*W)^2);
BOT=2.*W*sqrt(1+(W/2.)^2);
XMAG=TOP/BOT;
GAIN=20*log10(XMAG);
PHASE=57.3*atan2(W,29.1)-90.-57.3*atan2(W,2.);
count=count+1;
ArrayW(count)=W;
ArrayGAIN(count)=GAIN;
ArrayPHASE(count)=PHASE;
end
figure
semilogx(ArrayW,ArrayGAIN),grid
xlabel('Frequency (Rad/Sec)')
ylabel('Gain (Db)')
figure
semilogx(ArrayW,ArrayPHASE),grid
xlabel('Frequency (Rad/Sec)')
ylabel('Phase (Deg)')
clc
output=[ArrayW',ArrayGAIN',ArrayPHASE'];
save datfil.txt output -ascii
disp 'simulation finished' 	
 	