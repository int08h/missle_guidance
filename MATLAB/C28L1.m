clear
count=0;
VC=4000.;
BETA=161;
VM=3000.;
XNP=3.;
TAU=.2;
RUN=1000;
XLIM=999999999.;
for TF=.2:.2:10,
	Z1=0.;
	for I=1:RUN
		SUM=rand(1);
		TSTART=TF*SUM;
		PZ=rand(1);
		PZ=PZ-.5;
		if PZ > 0
			COEF=1;
		else
			COEF=-1;
		end;
		Y=0.;
		YD=0.;
		T=0.;
		H=.01;
		S=0.;
		XNC=0.;
		XNL=0.;
		while T <= (TF - 1e-5)
			YOLD=Y;
			YDOLD=YD;
			XNLOLD=XNL;
			STEP=1;
			FLAG=0;
			while STEP <=1
				if FLAG==1
					Y=Y+H*YD;
					YD=YD+H*YDD;
					XNL=XNL+H*XNLD;
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
				if T<TSTART
					XNT=0.;
				else
					XNT=COEF*BETA;
				end
				YDD=XNT-XNL;
				FLAG=1;
			end;
			FLAG=0;
			Y=.5*(YOLD+Y+H*YD);
			YD=.5*(YDOLD+YD+H*YDD);
			XNL=.5*(XNLOLD+XNL+H*XNLD);
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

