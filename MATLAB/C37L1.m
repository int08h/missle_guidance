clear
n=0;
IGUID=0;
ICHOICE=0;
GAMDEG=30.;
GAMIC=GAMDEG/57.3;
DELT=30.;
TBEG=0.;
TEND=TBEG+DELT;
XNP=3.;
RM1IC=0.;
RM2IC=0.;
RT1IC=10000.*3.28;
RT2IC=0.;
VM=250.*3.28;
XNCLIMG=10.;
GAMFDEG=-90.;
H=.0001;
XNCLIM=32.2*XNCLIMG;
RM1=RM1IC;
RM2=RM2IC;
RT1=RT1IC;
RT2=RT2IC;
VT1=0.;
VT2=0.;
T=0.;
S=0.;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1^2+RTM2^2);
XLAM=atan2(RTM2,RTM1);
VM1=VM*cos(GAMIC);
VM2=VM*sin(GAMIC);
VTM1=VT1-VM1;
VTM2=VT2-VM2;
VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
GAMF=GAMFDEG/57.3;
BIASDEG=(-GAMFDEG*(XNP-1.)+XNP*XLAM*57.3-GAMDEG)/DELT;
BIAS=BIASDEG/57.3;
X=0.;
while VC >= 0 
    if RTM < 1000
        H=.00001;
    else
        H=.0001;
    end
    RM1OLD=RM1;
    RM2OLD=RM2;
    VM1OLD=VM1;
    VM2OLD=VM2;
    XOLD=X;
    STEP=1;
    FLAG=0;
    while STEP <=1
        if FLAG==1
                STEP=2;
                RM1=RM1+H*VM1;
                RM2=RM2+H*VM2;
                VM1=VM1+H*AM1;
                VM2=VM2+H*AM2;
                X=X+H*XD;       
                T=T+H;
        end
        GAM=atan2(VM2,VM1);
        VM=sqrt(VM1^2+VM2^2);
        RTM1=RT1-RM1;
        RTM2=RT2-RM2;
        RTM=sqrt(RTM1^2+RTM2^2);
        VTM1=VT1-VM1;
        VTM2=VT2-VM2;
        VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
        XLAM=atan2(RTM2,RTM1);
        XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
        if IGUID==0
            TGO=RTM/VC;
            XNC=4.*VC*XLAMD+2.*VC*(XLAM-GAMF)/TGO;
        else
            if T<TBEG
                GAMD=XNP*XLAMD;
                XNC=VM*GAMD;
            elseif T<TEND
                GAMD=XNP*XLAMD+BIAS;
                XNC=VM*GAMD;
            else
                GAMD=XNP*XLAMD;
                XNC=VM*GAMD;
            end
        end
        if XNC>XNCLIM
            XNC=XNCLIM;
        end
        if XNC<-XNCLIM
            XNC=-XNCLIM;
        end
        if ICHOICE==0
            AM1=-XNC*sin(XLAM);
            AM2=XNC*cos(XLAM);
        else
            AM1=-XNC*sin(GAM);
            AM2=XNC*cos(GAM);
        end
        XD=XNC*XNC;
        FLAG=1;
    end
    FLAG=0;
    RM1=.5*(RM1OLD+RM1+H*VM1);
    RM2=.5*(RM2OLD+RM2+H*VM2);
    VM1=.5*(VM1OLD+VM1+H*AM1);
    VM2=.5*(VM2OLD+VM2+H*AM2);
    X=.5*(XOLD+X+H*XD);
    S=S+H;
    if S>=.09999    
        S=0.;
        n=n+1;
        RT1K=RT1/3280.;
        RT2K=RT2/3280.;
        RM1K=RM1/3280.;
        RM2K=RM2/3280.;
        XNCG=XNC/32.2;
        GAMDEG=GAM*57.3;
        VMM=VM/3.28;
        XM=X/(3.28*3.28);
        ArrayT(n)=T;
        ArrayRT1K(n)=RT1K;
        ArrayRT2K(n)=RT2K;
        ArrayRM1K(n)=RM1K;
        ArrayRM2K(n)=RM2K;
        ArrayXNCG(n)=XNCG;
        ArrayGAMDEG(n)=GAMDEG;
        ArrayXM(n)=XM;
    end
end
figure
plot(ArrayRT1K,ArrayRT2K,ArrayRM1K,ArrayRM2K),grid
title('Engagement Geometry')
xlabel('Downrange (Kft) ')
ylabel('Altitude (Kft)')
figure
plot(ArrayT,ArrayXNCG),grid
title('Commanded Acceleration')
xlabel('Time (Sec) ')
ylabel('XNC (G)')
figure
plot(ArrayT,ArrayGAMDEG),grid
title('Flight Path Angle')
xlabel('Time (Sec) ')
ylabel('GAM (Deg)')
clc
output=[ArrayT',ArrayRT1K',ArrayRT2K',ArrayRM1K',ArrayRM2K',...
		ArrayXNCG',ArrayGAMDEG'];
save datfil.txt output -ascii
disp '*** Simulation Complete'
RTM=sqrt(RTM1^2+RTM2^2)
VM=sqrt(VM1^2+VM2^2)
