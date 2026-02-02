count=0;
QACT=1;
QACT=0;
WACT=100.;
ZACT=.65;
XKD=9000.;
WRR=2.;
C1=3.;
C2=.103;
T=0;
H=.0002;
S=0.;
PHI=10.;
PHID=0.;
DEL=0.;
DELD=0.
while ~(T>.5)
	S=S+H;
	PHIOLD=PHI;
	PHIDOLD=PHID;
	DELOLD=DEL;
	DELDOLD=DELD;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			STEP=2;	
			PHI=PHI+H*PHID;
			PHID=PHID+H*PHIDD;
			DEL=DEL+H*DELD;
			DELD=DELD+H*DELDD;
			T=T+H;
		end
		DELC=-C1*PHI-C2*PHID;
		DELDD=WACT*WACT*(DELC-DEL-2.*ZACT*DELD/WACT);
        if QACT==1
            E=DEL;
        else
            E=DELC
        end
		PHIDD=XKD*E-WRR*PHID;
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
