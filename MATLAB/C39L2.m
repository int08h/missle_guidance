clear
VC=4000.;
GAM=.0001;
APN=2;
XNT=161;
YIC=0.;
VM=3000.;
HEDEG=0.;
TAU=.5;
WZ=10;
XNP=3.;
n=0.;
for TF=.1:.1:10
    Y=0;
    YD=0;
    E=0.;
    T=0.;
    H=.01;
    XNL=0;
    XNCLIM=9999999;
    while T<=(TF-1e-5)
        YOLD=Y;
        YDOLD=YD;
        EOLD=E;
        STEP=1;
        FLAG=0;
        while STEP<=1
      		if FLAG==1
         		STEP=2;
         		Y=Y+H*YD; 	
         		YD=YD+H*YDD;
         		E=E+H*ED;
         		T=T+H;
      		end
      		TGO=TF-T+.00001;
            XLAMD=(Y+YD*TGO)/(VC*TGO*TGO);
            if APN==0
                XNC=XNP*(Y+YD*TGO)/TGO^2;
            elseif APN==1
                 X=TGO/TAU;
                TOP=6.*X*X*(exp(-X)-1.+X);
                BOT1=2*X*X*X+3.+6.*X-6.*X*X;
                BOT2=-12.*X*exp(-X)-3*exp(-2.*X);
                XNPP=TOP/(.0001+BOT1+BOT2);
                XNEW=XNPP*XNL*(exp(-X)+X-1.)/(X*X);
                XNC=XNPP*VC*XLAMD+.5*XNPP*XNT-XNEW;
            else
                 XS=TGO/TAU;
                TEMP1=TGO*TGO*TAU*(exp(-XS)-1.+XS);
                TOP=-(TGO^3)/(TAU*WZ)+(1.+1./(TAU*WZ))*TEMP1;
                TEMP2=.5*(1.-3./(TAU*WZ))+XS*(1.+1./(TAU*WZ))-XS*XS;
                TEMP3=-2.*XS*exp(-XS);
                TEMP4=2.*exp(-XS)/(TAU*WZ)-.5*exp(-2.*XS)*(1.+1./(TAU*WZ));
	     BOT=GAM+TGO*TGO*TGO/3.+(1.+1./(TAU*WZ))*TAU^3*(TEMP2+TEMP3+TEMP4);            
	        XNPP=TOP/BOT;
	        C1TH=XNPP/(TGO*TGO);
	        C2TH=XNPP/TGO;
	        C3TH=.5*XNPP;
	        C4TH=-XNPP*(exp(-XS)+XS-1.)/(XS*XS);
	        XNC=C1TH*Y+C2TH*YD+C3TH*XNT+C4TH*E;
            end
      		if XNC>XNCLIM
         		XNC=XNCLIM;
      		end
      		if XNC<-XNCLIM
         		XNC=-XNCLIM;
      		end
            ED=(XNC-E)/TAU;
            XNL=E-ED/WZ;
           YDD=XNT-XNL;
      		FLAG=1;
        end
        FLAG=0;
        Y=.5*(YOLD+Y+H*YD);
        YD=.5*(YDOLD+YD+H*YDD);
        E=.5*(EOLD+E+H*ED);
    end
    n=n+1;
    ArrayTF(n)=TF;
    ArrayY(n)=Y;
end
figure
plot(ArrayTF,ArrayY),grid
xlabel('Flight Time (Sec)')
ylabel('Miss (ft)')
clc
output=[ArrayTF',ArrayY',];
save datfil.txt output  -ascii
disp 'simulation finished'
