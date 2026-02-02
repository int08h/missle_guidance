clear
n=0;
ZACT=.7;
WACT=150;
K3=-1.89;
TA=.457;
ZAF=.058;
WAF=25.3;
KR=.1;
PI=3.1416;
H=.0001;
for I=2:160
	W=10^(.025*I-1);
	PERIOD=2.*PI/W;
	T=0.;
	S=0.;
	E=0.;
	ED=0.;
	DEL=0.;
	DELD=0.;
	P=0.;
	Q=0;
	PPREV=0;
	QPREV=0;
	DELP=0;
	DELQ=0;
	DELPOLD=0;
	DELQOLD=0;
	DELDELP=100;
	DELDELQ=100; 
	while ~(T>20. & abs(DELDELP)<.0001)
		EOLD=E;
		EDOLD=ED;
		DELOLD=DEL;
		DELDOLD=DELD;
		POLD=P;
		QOLD=Q;
		STEP=1;
		FLAG=0;
		while STEP<=1
      		if FLAG==1
         		STEP=2;
				E=E+H*ED;
				ED=ED+H*EDD;
				DEL=DEL+H*DELD;
				DELD=DELD+H*DELDD;
				P=P+H*PD;
				Q=Q+H*QD;
				T=T+H;
			end
			X=-sin(W*T);
			DELDD=WACT*WACT*(X-DEL-2.*ZACT*DELD/WACT);
			EDD=WAF*WAF*(DEL-E-2.*ZAF*ED/WAF);
			Y=KR*K3*(E+TA*ED);
			PD=Y*sin(W*T);
			QD=Y*cos(W*T);
			FLAG=1;
   		end
   		FLAG=0; 
		E=.5*(EOLD+E+H*ED);
		ED=.5*(EDOLD+ED+H*EDD);
		DEL=.5*(DELOLD+DEL+H*DELD);
		DELD=.5*(DELDOLD+DELD+H*DELDD);
		P=.5*(POLD+P+H*PD);
		Q=.5*(QOLD+Q+H*QD);
		S=S+H;
		if (S>=(PERIOD-.0001))
		S=0.;
			DELP=P-PPREV;
			DELQ=Q-QPREV;
			PPREV=P;
			QPREV=Q;
			DELDELP=DELPOLD-DELP;
			DELDELQ=DELQOLD-DELQ;
			DELPOLD=DELP;
			DELQOLD=DELQ;
		end
	end
	PHASE=57.3*atan2(DELQ,DELP);
	if PHASE>90.
		PHASE=PHASE-360;
	end
	GAIN=10.*log10((DELP^2+DELQ^2)*W*W/(PI*PI));
	n=n+1;
	ArrayW(n)=W;
    ArrayPHASE(n)=PHASE;
    ArrayGAIN(n)=GAIN; 
end
figure
plot(ArrayW,ArrayGAIN),grid
xlabel('Frequency (r/s)')
ylabel('Gain (db)')
figure
plot(ArrayW,ArrayPHASE),grid
xlabel('Frequency (r/s)')
ylabel('Phase (deg)')
clc
output=[ArrayW',ArrayGAIN',ArrayPHASE'];
save datfil output  -ascii
disp 'simulation finished'
