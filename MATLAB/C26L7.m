clear
count=0;
for I=2:160
	W=10^(.025*I-1);
	TOP=13500*sqrt(1+(W/29.1)^2);
	BOT=W*sqrt(1+(W/2.)^2);
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
axis([.1 1000 -60 40])
figure
semilogx(ArrayW,ArrayPHASE),grid
xlabel('Frequency (Rad/Sec)')
ylabel('Phase (Deg)')
axis([.1 1000 -400 100])
clc
output=[ArrayW',ArrayGAIN',ArrayPHASE'];
save datfil.txt output -ascii
disp 'simulation finished' 	