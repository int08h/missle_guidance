clear
count=0;
QPN=1;
TAU=1.;
W=3.;
AT=193.2;
VT=1000.;
VM=3000.;
XNP=3.;
XNCLIM=9999999999.;
for RT3IC=40000:-500:500
	RM1=0.;
	RM2=10000.;
	RM3=0.;
	RT1=0.;
	RT2=10000.;
	RT3=RT3IC;
	VT1=-AT/W;
	VT2=0.;
	VT3=-VT;
	T=0.;
	S=0.;
	RTM1=RT1-RM1;
	RTM2=RT2-RM2;
	RTM3=RT3-RM3;
	RTM=sqrt(RTM1^2+RTM2^2+RTM3^2);
	VM1=0.;
	VM2=0.;
	VM3=VM;
	VTM1=VT1-VM1;
	VTM2=VT2-VM2;
	VTM3=VT3-VM3;
	VC=-(RTM1*VTM1+RTM2*VTM2+RTM3*VTM3)/RTM;
	AM1=0.;
	AM2=0.;
	AM3=0.;
 	while VC >= 0 
 		if RTM < 1000
      			H=.0002;
   		else
      			H=.01;
   		end 
		RT1OLD=RT1;
		RT2OLD=RT2;
		RT3OLD=RT3;
		RM1OLD=RM1;
		RM2OLD=RM2;
		RM3OLD=RM3;
		VM1OLD=VM1;
		VM2OLD=VM2;
		VM3OLD=VM3;
		VT1OLD=VT1;
		VT2OLD=VT2;
		VT3OLD=VT3;
		AM1OLD=AM1;
		AM2OLD=AM2;
		AM3OLD=AM3;
		STEP=1;
		FLAG=0;
		while STEP<=1
			if FLAG==1 
				STEP=2;
 				RT1=RT1+H*VT1;
				RT2=RT2+H*VT2;
				RT3=RT3+H*VT3;
				RM1=RM1+H*VM1;
				RM2=RM2+H*VM2;
				RM3=RM3+H*VM3;
				VM1=VM1+H*AM1;
				VM2=VM2+H*AM2;
				VM3=VM3+H*AM3;
				VT1=VT1+H*AT1;
				VT2=VT2+H*AT2;
				VT3=VT3+H*AT3;
				AM1=AM1+H*AM1D;
				AM2=AM2+H*AM2D;
				AM3=AM3+H*AM3D;
				T=T+H;
			end
			RTM1=RT1-RM1;
			RTM2=RT2-RM2;
			RTM3=RT3-RM3;
			RTM=sqrt(RTM1^2+RTM2^2+RTM3^2);
			VTM1=VT1-VM1;
			VTM2=VT2-VM2;
			VTM3=VT3-VM3;
			VC=-(RTM1*VTM1+RTM2*VTM2+RTM3*VTM3)/RTM;
			TGO=RTM/VC;
			AT1=AT*sin(W*T);
			AT2=AT*cos(W*T);
			AT3=0.;
			if QPN==1 
				ZEM1=RTM1+VTM1*TGO;
				ZEM2=RTM2+VTM2*TGO;
				ZEM3=RTM3+VTM3*TGO;
			else
				AT1D=AT*W*cos(W*T);
				AT2D=-AT*W*sin(W*T);
				AT3D=0.;
				X=TGO/TAU;
				TOP=6.*X*X*(exp(-X)-1.+X);
				BOT1=2*X*X*X+3.+6.*X-6.*X*X;
				BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
				XNP=TOP/(.0001+BOT1+BOT2);
				ZEM1=RTM1+VTM1*TGO+AT1*(1.-cos(W*TGO))/W^2+...
					(AT1D*(W*TGO-sin(W*TGO))/W^3)-...
					AM1*TAU*TAU*(exp(-X)+X-1.);
				ZEM2=RTM2+VTM2*TGO+AT2*(1.-cos(W*TGO))/W^2+...
					(AT2D*(W*TGO-sin(W*TGO))/W^3)-...
					AM2*TAU*TAU*(exp(-X)+X-1.);
				ZEM3=RTM3+VTM3*TGO+AT3*(1.-cos(W*TGO))/W^2+...
					(AT3D*(W*TGO-sin(W*TGO))/W^3)-...
					AM3*TAU*TAU*(exp(-X)+X-1.);
			end
			ZEMDOTRTM=(ZEM1*RTM1+ZEM2*RTM2+ZEM3*RTM3)/RTM;
			ZEMPER1=ZEM1-ZEMDOTRTM*RTM1/RTM;
			ZEMPER2=ZEM2-ZEMDOTRTM*RTM2/RTM;
			ZEMPER3=ZEM3-ZEMDOTRTM*RTM3/RTM;
			AM1P=XNP*ZEMPER1/(TGO*TGO);
			AM2P=XNP*ZEMPER2/(TGO*TGO);
			AM3P=XNP*ZEMPER3/(TGO*TGO);
			AM1D=(AM1P-AM1)/TAU;
			AM2D=(AM2P-AM2)/TAU;
			AM3D=(AM3P-AM3)/TAU;
			if AM1>XNCLIM
				AM1=XNCLIM;
			end
			if AM1<-XNCLIM 
				AM1=-XNCLIM;
			end
			if AM2>XNCLIM 
				AM2=XNCLIM;
			end
			if AM2<-XNCLIM
				AM2=-XNCLIM;
			end
			if AM3>XNCLIM
				AM3=XNCLIM;
			end
			if AM3<-XNCLIM 
				AM3=-XNCLIM;
			end
			XNCG=sqrt(AM1^2+AM2^2+AM3^2)/32.2;
			FLAG=1;
		end
		FLAG=0;
 		RT1=.5*(RT1OLD+RT1+H*VT1);
		RT2=.5*(RT2OLD+RT2+H*VT2);
		RT3=.5*(RT3OLD+RT3+H*VT3);
		RM1=.5*(RM1OLD+RM1+H*VM1);
		RM2=.5*(RM2OLD+RM2+H*VM2);
		RM3=.5*(RM3OLD+RM3+H*VM3);
		VM1=.5*(VM1OLD+VM1+H*AM1);
		VM2=.5*(VM2OLD+VM2+H*AM2);
		VM3=.5*(VM3OLD+VM3+H*AM3);
		VT1=.5*(VT1OLD+VT1+H*AT1);
		VT2=.5*(VT2OLD+VT2+H*AT2);
		VT3=.5*(VT3OLD+VT3+H*AT3);
		AM1=.5*(AM1OLD+AM1+H*AM1D);
		AM2=.5*(AM2OLD+AM2+H*AM2D);
		AM3=.5*(AM3OLD+AM3+H*AM3D);
		S=S+H;
		if S>=.09999
			S=0.;
			RT1K=RT1/1000.;
			RT2K=RT2/1000.;
			RT3K=RT3/1000.;
			RM1K=RM1/1000.;
			RM2K=RM2/1000.;
			RM3K=RM3/1000.;

		end
	end
	count=count+1;
	ArrayT(count)=T;
	ArrayRTM1(count)=RTM1;
	ArrayRTM2(count)=RTM2;
	ArrayRTM3(count)=RTM3;
	ArrayRTM(count)=RTM;
end
output=[ArrayT',ArrayRTM1',ArrayRTM2',ArrayRTM3',ArrayRTM'];
save datfil.txt output -ascii
disp 'simulation finished'
clc
figure
plot(ArrayT,ArrayRTM1),grid
xlabel('Time (s) ')
ylabel('RTM1 (ft)')	
axis([0 10 -30 30])
figure
plot(ArrayT,ArrayRTM2),grid
xlabel('Time (s) ')
ylabel('RTM2 (ft)')	
axis([0 10 -30 30])
figure
plot(ArrayT,ArrayRTM3),grid
xlabel('Time (s) ')
ylabel('RTM3 (ft)')	
axis([0 10 -30 30])
figure
plot(ArrayT,ArrayRTM),grid
xlabel('Time (s) ')
ylabel('RTM (ft)')	
axis([0 10 0 30])

