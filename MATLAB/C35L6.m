clear
count=0;
TAU=1.;
W=20.;
Z=.7;
WZ=5.;
APN=3;
XNT=16.1;
TS=.01;
GAM=.00001;
XLIM=9999999.;
TFMAX=10.;
XNP=3.;
VC=4000.;
VM=3000.;
if APN==3
	[C1,C2,C3,C4,C5,C6]=GENERATEGAINS(TAU,W,Z,WZ,GAM,TFMAX,TS);
end
for TF=.1:.1:10.
	E=0.;
	ED=0.;
	EDD=0.;
	T=0;
	H=.0001;
	S=0.;
	Y=0.;
	YD=0.;
	XNC=0.;
	RTM=VC*TF;
	while T< (TF-.00001)
		S=S+H;
		EOLD=E;
		EDOLD=ED;
		EDDOLD=EDD;
		YOLD=Y;
		YDOLD=YD;
		STEP=1;
		FLAG=0;
		while STEP<=1
			if FLAG==1
         			STEP=2	;
         			E=E+H*ED;
				ED=ED+H*EDD;
				EDD=EDD+H*EDDD;
				Y=Y+H*YD;
				YD=YD+H*YDD;
				T=T+H;
			end
			TGO=TF-T+.0001;
			RTM=VC*TGO;
			XLAM=Y/RTM;
			XNCG=XNC/32.2;
 			EDDD=W*W*(XNC-E-(2.*Z/W+TAU)*ED-(2.*Z*TAU/W+1./W^2)...
					  *EDD)/TAU;
 			XNL=E-EDD/WZ^2;
			YDD=XNT-XNL;
			FLAG=1;
		end
		FLAG=0;
		E=.5*(EOLD+E+H*ED);
		ED=.5*(EDOLD+ED+H*EDD);
		EDD=.5*(EDDOLD+EDD+H*EDDD);
		Y=.5*(YOLD+Y+H*YD);
		YD=.5*(YDOLD+YD+H*YDD);
		if S>=(TS-.0001)
			S=0.;
			if APN==0
				XNC=XNP*(Y+YD*TGO)/(TGO*TGO);
				XNPP=XNP;
			elseif APN==1
				XNC=XNP*(Y+YD*TGO+.5*XNT*TGO*TGO)/(TGO*TGO);
				XNPP=XNP;
			elseif APN==2
				XS=TGO/TAU;
				TOP=6.*XS*XS*(exp(-XS)-1.+XS);
				BOT1=2*XS*XS*XS+3.+6.*XS-6.*XS*XS;
				BOT2=-12.*XS*exp(-XS)-3.*exp(-2.*XS);
				XNPP=TOP/(.0001+BOT1+BOT2);
				C1P=XNPP/(TGO*TGO);
				C2P=XNPP/TGO;
				C3P=.5*XNPP;
				C4P=-XNPP*(exp(-XS)+XS-1.)/(XS*XS);
				XNC=C1P*Y+C2P*YD+C3P*XNT+C4P*XNL;
			else
				JJ=fix(TGO/TS)+1;
				XNC=C1(JJ)*Y+C2(JJ)*YD+C3(JJ)*XNT+C4(JJ)*E...
					+C5(JJ)*ED+C6(JJ)*EDD;
     				XNPP=C2(JJ)*TGO;
     			end
			if XNC>XLIM
				XNC=XLIM;
			elseif XNC<-XLIM
				XNC=-XLIM;
			end
			XNCG=XNC/32.2;
		end
	end
	count=count+1;
	ArrayTF(count)=TF;
	ArrayY(count)=Y;
end
figure
plot(ArrayTF,ArrayY),grid
xlabel('Flight Time (s)')
ylabel('Miss (ft) ')
clc
output=[ArrayTF',ArrayY'];
save datfil.txt output -ascii
	

	
	
