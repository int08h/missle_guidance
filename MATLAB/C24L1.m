clear
n=0;
VM=3000.;
DEL=5./57.3;
ALT=0.;
A=1000.;
DIAM=1.;
FR=3.;
XL=20.;
CTW=0.;
CRW=6.;
HW=2.;
CTT=0.;
CRT=2.;
HT=2.;
XN=4.;
XCG=10.;
XHL=19.5;
WGT=1000.;
if ALT<=30000.
	RHO=.002378*exp(-ALT/30000.);
else
	RHO=.0034*exp(-ALT/22000.);
end
SWING=.5*HW*(CTW+CRW);
STAIL=.5*HT*(CTT+CRT);
SREF=3.1416*DIAM*DIAM/4.;
XLP=FR*DIAM;
SPLAN=(XL-XLP)*DIAM+1.33*XLP*DIAM/2.;
XCPN=2*XLP/3;
AN=.67*XLP*DIAM;
AB=(XL-XLP)*DIAM;
XCPB=(.67*AN*XLP+AB*(XLP+.5*(XL-XLP)))/(AN+AB);
XCPW=XLP+XN+.7*CRW-.2*CTW;
XMACH=VM/A;
XIYY=WGT*(3*((DIAM/2)^2)+XL*XL)/(12*32.2);
TMP1=(XCG-XCPW)/DIAM;
TMP2=(XCG-XHL)/DIAM;
TMP3=(XCG-XCPB)/DIAM;
TMP4=(XCG-XCPN)/DIAM;
B=sqrt(XMACH^2-1);
Q=.5*RHO*VM*VM;
THD=0;
ALF=0;
T=0;
H=.0025;
S=0.;
while T<1.99999
	THDOLD=THD;
	ALFOLD=ALF;
	STEP=1;
	FLAG=0;
   	while STEP<=1
		if FLAG==1
			STEP=2;
			THD=THD+H*THDD;
			ALF=ALF+H*ALFD;
			T=T+H;
		end
		CN=2*ALF+1.5*SPLAN*ALF*ALF/SREF+8*SWING*ALF/(B*SREF)+...
         		8*STAIL*(ALF+DEL)/(B*SREF);
		CM=2*ALF*TMP4+1.5*SPLAN*ALF*ALF*TMP3/SREF+...
			8*SWING*ALF*TMP1/(B*SREF)...
         		+8*STAIL*(ALF+DEL)*TMP2/(B*SREF);
		THDD=Q*SREF*DIAM*CM/XIYY;
		XNL=32.2*Q*SREF*CN/WGT;
		ALFD=THD-XNL/VM;
      		FLAG=1;
   	end
   	FLAG=0;
	THD=.5*(THDOLD+THD+H*THDD);
	ALF=.5*(ALFOLD+ALF+H*ALFD);
	S=S+H;
	if S>=.0099999
		S=0.;
		n=n+1;
		ArrayT(n)=T;
		ArrayXNLG(n)=XNL/32.2;
		ArrayALFDEG(n)=ALF*57.3;
	end
end
figure
plot(ArrayT,ArrayXNLG),grid
xlabel('Time (Sec)')
ylabel('Missile Acceleration (G)')
figure
plot(ArrayT,ArrayALFDEG),grid
xlabel('Time (Sec)')
ylabel('Angle of Attack (Deg)')
clc
output=[ArrayT',ArrayXNLG',ArrayALFDEG'];
save datfil.txt output  -ascii
disp 'simulation finished'
