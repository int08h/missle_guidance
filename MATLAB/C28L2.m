clear
count=0;
VC=4000.;
BETA=161;
VM=3000.;
XNP=3.;
TAU=.2;
RUN=1000;
XLIM=999999999.;
H=.01;
for TF=.2:.2:10,
	Z1=0.;
	for I=1:RUN
		PHI=BETA*BETA/TF;
		SIG=sqrt(PHI/H);
		Y=0.;
		YD=0.;
		T=0.;
		S=0.;
		XNC=0.;
		XNL=0.;
		XNT=0.;
		while T <= (TF - 1e-5)
			X=SIG*randn;
			YOLD=Y;
			YDOLD=YD;
			XNLOLD=XNL;
			XNTOLD=XNT;
			STEP=1;
			FLAG=0;
			while STEP <=1
				if FLAG==1
					Y=Y+H*YD;
					YD=YD+H*YDD;
					XNL=XNL+H*XNLD;
					XNT=XNT+H*XNTD;
					T=T+H;
					STEP=2;
				end;
				TGO=TF-T+.00001;
				RTM=VC*TGO;
				XLAM=Y/RTM;
				XLAMD=(RTM*YD+Y*VC)/(RTM^2);
				XNC=XNP*VC*XLAMD;
				if XNC>XLIM
					XNC=XLIM;
				elseif XNC<-XLIM
					XNC=-XLIM;
				end;
				XNLD=(XNC-XNL)/TAU;
				YDD=XNT-XNL;
				XNTD=X;
				FLAG=1;
			end;
			FLAG=0;
			Y=.5*(YOLD+Y+H*YD);
			YD=.5*(YDOLD+YD+H*YDD);
			XNL=.5*(XNLOLD+XNL+H*XNLD);
			XNT=.5*(XNTOLD+XNT+H*XNTD);
			S=S+H;
		end
		Z(I)=Y;
		Z1=Z(I)+Z1;
		XMEAN=Z1/I;
	end
	SIGMA=0.;
	Z1=0.;
	Z2=0.;
	for I=1:RUN,
		Z1=(Z(I)-XMEAN)^2+Z1;
		Z2=Z(I)^2+Z2;
		if I==1,
			SIGMA=0;
			RMS=0.;
		else
			SIGMA=sqrt(Z1/(I-1));
			RMS=sqrt(Z2/(I-1));
		end; 
    end
	count=count+1;
	ArrayTF(count)=TF;
	ArrayRMS(count)=RMS;
end;
figure
plot(ArrayTF',ArrayRMS'),grid
title('RMS miss for various flight times')
xlabel('Flight Time (S)')
ylabel('RMS MISS (Ft) ')
clc
output=[ArrayTF',ArrayRMS'];
save datfil.txt output -ascii
disp('Simulation Complete')

