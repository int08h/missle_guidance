clear
count=0;
ZACT=.7;
WACT=150.;
K3=-1.89;
TA=.457;
ZAF=.058;
WAF=25.3;
KR=.1;
for I=2:160
	W=10^(.025*I-1);
	XMAG1=sqrt(1+(W*TA)^2);
	XMAG2=sqrt((1-(W/WAF)^2)^2+(2*ZAF*W/WAF)^2);
	XMAG3=sqrt((1-(W/WACT)^2)^2+(2*ZACT*W/WACT)^2);
	GAIN=20*log10(-K3*KR*XMAG1/(XMAG2*XMAG3));
	PHASE1=57.3*atan2(W*TA,1.);
	PHASE2=57.3*atan2(2*ZAF*W/WAF,1-(W/WAF)^2);
	PHASE3=57.3*atan2(2*ZACT*W/WACT,1-(W/WACT)^2);
	PHASE=PHASE1-PHASE2-PHASE3;
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


