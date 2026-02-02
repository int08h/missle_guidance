clear
count=0;
XNT=161,;
XNP=3.;
TAU=.5;
TF=10.;
VM=2000.;
HEDEG=0.;
P=3.;
XL=P/2.;
TR=1.;
PI=3.1416;
ALF=PI*TR/(2.*XL);
W=2.*PI/P;
VC=3000.;
T=0.;
	S=0.
TP=T+.00001;
X1=0.;
X2=0;
X3=1.;
X4=0;
X5=0.;
X6=0.;
X7=0.;
X8=0.;
H=.01;
HE=HEDEG/57.3;
while TP<=(TF-1e-5)
    X1OLD=X1;
    X2OLD=X2;
    X3OLD=X3;
    X4OLD=X4;
    X5OLD=X5;
    X6OLD=X6;
    X7OLD=X7;
    X8OLD=X8;
    STEP=1;
    FLAG=0;
    while STEP<=1
        if FLAG==1
            STEP=2;
            X1=X1+H*X1D;
            X2=X2+H*X2D;
            X3=X3+H*X3D;
            X4=X4+H*X4D;
            X5=X5+H*X5D;
            X6=X6+H*X6D;
            X7=X7+H*X7D;
            X8=X8+H*X8D;
            TP=TP+H;
        end
        TGO=TP+.00001;
        X1D=X5;
        X2D=X3;
        Y1=X4-XNP*VC*X2;
        X3D=Y1/(VC*TGO*TAU);
        X4D=-Y1/TAU;
        X5D=X2-W*W*X1;
        X6D=X7;
        X7D=X2-9.*W*W*X6;
        Y2=X1*sin(ALF)+X6*sin(3.*ALF)/9.;
        PZ=4.*W*Y2/(PI*ALF);
        X8D=PZ^2;
        FLAG=1;
    end
    FLAG=0;
    X1=(X1OLD+X1)/2+.5*H*X1D;
    X2=(X2OLD+X2)/2+.5*H*X2D;
    X3=(X3OLD+X3)/2+.5*H*X3D;
    X4=(X4OLD+X4)/2+.5*H*X4D;
    X5=(X5OLD+X5)/2+.5*H*X5D;
    X6=(X6OLD+X6)/2+.5*H*X6D;
    X7=(X7OLD+X7)/2+.5*H*X7D;
    X8=(X8OLD+X8)/2+.5*H*X8D;
    S=S+H;
    if S>=.0999
        S=0.;
        count=count+1;
        XMNT=4.*W*XNT*Y2/(PI*ALF);
        XMUDNT=XNT*sqrt(X8/TGO);
        ArrayTP(count)=TP;
        ArrayXMNT(count)=XMNT;
        ArrayXMUDNT(count)=XMUDNT;
    end
end
figure
plot(ArrayTP,ArrayXMNT),grid
xlabel('Flight Time (Sec)')
ylabel('Target Maneuver Miss (Ft)')
figure
plot(ArrayTP,ArrayXMUDNT),grid
xlabel('Flight Time (Sec)')
ylabel('Random Target Maneuver Miss (Ft)')
clc
output=[ArrayTP',ArrayXMNT',ArrayXMUDNT'];
save datfil output  -ascii
disp 'simulation finished'
