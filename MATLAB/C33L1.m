clear
count=0;
QPN=1;
IMVR=1;
VC=4000.;
Y=0.;
VM=3000.;
XNTIC=96.6;
HEDEG=0.;
W=3.;
TF=10.;
XNP=3.;
XNPP=10.;
XNCDLIM=999999.;
YD=-VM*HEDEG/57.3;
T=0.;
H=.001;
S=0.;
XNC=XNP*YD/TF;
J=0.;
while T<=(TF - 1e-5)
    YOLD=Y;
    YDOLD=YD;
    XNCOLD=XNC;
    JOLD=J;
    STEP=1;
    FLAG=0;
    while STEP <=1
        if FLAG==1
            Y=Y+H*YD;
            YD=YD+H*YDD;
            XNC=XNC+H*XNCD;
            J=J+H*JD;
            T=T+H;
            STEP=2;
        end
        TGO=TF-T+.00001;
        if IMVR==1
            XNT=XNTIC;
        else
            if sin(W*T)>0.
                XNT=XNTIC;
            end
            if sin(W*T)<=0
                XNT=-XNTIC;
            end
        end
        YDD=XNT-XNC;
        if QPN==1
            XNCD=2.*XNP*(Y+YD*TGO+.5*YDD*TGO*TGO)/TGO^3;
        else
            XNCD=XNPP*(Y+YD*TGO+.5*YDD*TGO*TGO)/TGO^3;
        end
        if XNCD>XNCDLIM
            XNCD=XNCDLIM;
        end
        if XNCD<-XNCDLIM
            XNCD=-XNCDLIM;
        end
        JD=XNCD^2;
        FLAG=1;
    end
    FLAG=0;
    Y=.5*(YOLD+Y+H*YD);
    YD=.5*(YDOLD+YD+H*YDD);
    XNC=.5*(XNCOLD+XNC+H*XNCD);
    J=.5*(JOLD+J+H*JD);
    S=S+H;
	if S>=.09999
        S=0.;
        count=count+1;
        ArrayT(count)=T;
        ArrayY(count)=Y;
        ArrayXNCG(count)=XNC/32.2;
        ArrayXNCDG(count)=XNCD/32.2;
        ArrayXNTG(count)=XNT/32.2;
	end
end
figure
plot(ArrayT,ArrayY),grid
xlabel('Time (s)')
ylabel('Y (ft)  ')
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (s)')
ylabel('XNC (g)  ')
figure
plot(ArrayT,ArrayXNCDG),grid
xlabel('Time (s)')
ylabel('XNCD (g/s)  ')
figure
plot(ArrayT,ArrayXNTG),grid
xlabel('Time (s)')
ylabel('XNT (g)  ')
clc
output=[ArrayT',ArrayY',ArrayXNCG',ArrayXNCDG',ArrayXNTG'];
save datfil.txt output -ascii
disp 'simulation finished'
Y
