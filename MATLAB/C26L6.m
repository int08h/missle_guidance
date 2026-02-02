clear
count=0;
WG=200.;
ZG=.5;
WACT=100.;
ZACT=.65;
XKD=9000.;
WR=2.;
C1=3.;
C2=.127;
C3=8.84;
C4=.031;
C5=.0152;
C6=.0000768;
T=0;
H=.0002;
S=0.;
PHI=10.;
PHID=0.;
DEL=0.;
DELD=0.;
PHIM=10.;
PHIMD=0.;
PHIMDD=0.;
while ~(T>.5)
	S=S+H;
	PHIOLD=PHI;
	PHIDOLD=PHID;
	DELOLD=DEL;
	DELDOLD=DELD;
	PHIMOLD=PHIM;
	PHIMDOLD=PHIMD;
	PHIMDDOLD=PHIMDD;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			STEP=2;	
			PHI=PHI+H*PHID;
			PHID=PHID+H*PHIDD;
			DEL=DEL+H*DELD;
			DELD=DELD+H*DELDD;
			PHIM=PHIM+H*PHIMD;
			PHIMD=PHIMD+H*PHIMDD;
			PHIMDD=PHIMDD+H*PHIMDDD;
			T=T+H;
		end
		DELC=-C1*PHI-C2*PHID-C3*DEL-C4*DELD-C5*PHIMD-C6*PHIMDD;
		DELDD=WACT*WACT*(DELC-DEL-2.*ZACT*DELD/WACT);
		E=DEL;
		PHIDD=XKD*E-WR*PHID;
		PHIMDDD=WG*WG*(PHID-PHIMD-2.*ZG*PHIMDD/WG);
		FLAG=1;
	end
	FLAG=0;
	PHI=.5*(PHIOLD+PHI+H*PHID);
	PHID=.5*(PHIDOLD+PHID+H*PHIDD);
	if S>=.0-9999	
		S=0.;
		count=count+1;
		ArrayT(count)=T;
		ArrayPHI(count)=PHI;
		ArrayPHID(count)=PHID;
		ArrayDEL(count)=DEL;
	end
end
figure
plot(ArrayT,ArrayPHI),grid
xlabel('Time (s) ')
ylabel('PHI (deg)')
clc
figure
plot(ArrayT,ArrayPHID),grid
xlabel('Time (s) ')
ylabel('PHID (deg/s)')
figure
plot(ArrayT,ArrayDEL),grid
xlabel('Time (s) ')
ylabel('DEL (deg)')
clc
clc
output=[ArrayT',ArrayPHI',ArrayPHID',ArrayDEL'];
save datfil.txt output -ascii
disp '*** Simulation Complete'
