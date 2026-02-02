count=0;
RUN=100;
TAU=.05;
ORDER=6;
TSTART=0.;
PZ=3.;
TR=.5;
XL=PZ/2.;
PI=3.14159;
W=2.*PI/PZ;
WH=W;
WREAL=W;
VC=2.5*3280.;
XNT=322.;
XNTREAL=XNT;
TS=.01;
VM=2.*3280.;
HEDEGFIL=20.;
XNP=3.;
SIGRIN=.0001;
SIGGL=0.;
RA=21000.;
SRN=0.;
TF=10.;
QPERFECT=0;
PHASE=0./57.3;
X=WH*TS;
Y=0;
YD=0.;
XLIM=1288.;

        TSTART=0.;
        XNTIC=XNT;
        XL=PZ/2.;
        X=PZ/2.-TR;
        PI=3.1416;
        ALF=PI*TR/(2.*XL);
        W=2.*PI/PZ;
        PHIS=XNT*XNT/6.;
        RTM=VC*TF;
        SIGNOISE=sqrt(SIGRIN^2+(SIGGL/RTM)^2+(SRN*RTM*RTM/(RA*RA))^2);
        SIGPOS=RTM*SIGNOISE;
        SIGN2=SIGPOS^2;
        PHI=zeros(ORDER);
        P=zeros(ORDER);
        QC=zeros(ORDER);
        Q=zeros(ORDER);
        F=zeros(ORDER);
        IDNP=eye(ORDER);
        F(1,2)=1.;
        F(2,3)=1.;
        F(2,5)=1.;
        F(3,4)=1.;
        F(4,3)=-W*W;
        F(5,6)=1.;
        F(6,5)=-9.*W*W;
        PHI=IDNP+F*TS+.5*F*F*TS*TS;
        TMPA=4.*W*sin(ALF)/(PI*ALF);
        TMPB=4.*W*sin(3*ALF)/(3.*PI*ALF);
        QC(4,4)=PHIS*TMPA*TMPA;
        QC(4,6)=PHIS*TMPA*TMPB;
        QC(6,4)=PHIS*TMPA*TMPB;
        QC(6,6)=PHIS*TMPB*TMPB;
        Q=PHI*QC*PHI'*TS;
        P(1,1)=SIGN2;
        P(2,2)=(VM*HEDEGFIL/57.3)^2;
        P(3,3)=322.^2;
        P(4,4)=(W*322.)^2;
        P(5,5)=(322.)^2;
        P(6,6)=(W*322.)^2;
        HMAT=[1 0 0 0 0 0];
        T=0.;
        H=.001;
        S=0.;
        XNC=0.;
        XNL=0.;
        XLAM=Y/RTM;
        YTDD=0;
        YTDDD=0;
        YH=0.;
        YDH=0.;
        X1H=0.;
        X2H=0.;
        X3H=0.;
        X4H=0.;
        Y=0.;
        YD=0.;
        while T<=(TF-.0001)
            YOLD=Y;
            YDOLD=YD;
            XNLOLD=XNL;
            STEP=1;
            FLAG=0;
            while STEP <=1
                if FLAG==1
                    STEP=2;
                    Y=Y+H*YD;
                    YD=YD+H*YDD;
                    XNL=XNL+H*XNLD;
                    T=T+H;
                end
                TGO=TF-T+.000001;
                RTM=VC*TGO;
                XLAM=Y/(VC*TGO);
                if T<TSTART
                    YTDD=0.;
                elseif (T-TSTART)<PZ
                    if (T-TSTART)<TR/2.
                        YTDD=2.*XNT*(T-TSTART)/TR;
                    elseif (T-TSTART)<(TR/2.+X)
                        YTDD=XNT;
                    elseif (T-TSTART)<(3.*TR/2.+X)
                        YTDD=-2.*XNT*(T-TSTART)/TR+2.*XNT+2.*XNT*X/TR;
                    elseif (T-TSTART)<(3.*TR/2.+2.*X)
                        YTDD=-XNT;
                    else
                        YTDD=2.*XNT*(T-TSTART)/TR-4.*XNT-4.*XNT*X/TR;
                    end
                elseif (T-TSTART)<2.*PZ
                    TSTAR=(T-TSTART)-PZ;
                    if TSTAR<TR/2.
                        YTDD=2.*XNT*TSTAR/TR;
                    elseif TSTAR<(TR/2.+X)
                        YTDD=XNT;
                    elseif TSTAR<(3.*TR/2.+X)
                        YTDD=-2.*XNT*TSTAR/TR+2.*XNT+2.*XNT*X/TR;
                    elseif TSTAR<(3.*TR/2.+2.*X)
                        YTDD=-XNT;
                    else
                        YTDD=2.*XNT*TSTAR/TR-4.*XNT-4.*XNT*X/TR;
                    end
                elseif (T-TSTART)<3.*PZ
                    TSTAR=(T-TSTART)-2.*PZ;
                    if TSTAR<TR/2.
                        YTDD=2.*XNT*TSTAR/TR;
                    elseif TSTAR<(TR/2.+X)
                        YTDD=XNT;
                    elseif TSTAR<(3.*TR/2.+X)
                        YTDD=-2.*XNT*TSTAR/TR+2.*XNT+2.*XNT*X/TR;
                    elseif TSTAR<(3.*TR/2.+2.*X)
                        YTDD=-XNT;
                    else
                        YTDD=2.*XNT*TSTAR/TR-4.*XNT-4.*XNT*X/TR;
                    end
                else
                    TSTAR=(T-TSTART)-3.*PZ;
                    if TSTAR<TR/2.
                        YTDD=2.*XNT*TSTAR/TR;
                    elseif TSTAR<(TR/2.+X)
                        YTDD=XNT;
                    elseif TSTAR<(3.*TR/2.+X)
                        YTDD=-2.*XNT*TSTAR/TR+2.*XNT+2.*XNT*X/TR;
                    elseif TSTAR<(3.*TR/2.+2.*X)
                        YTDD=-XNT;
                    else
                        YTDD=2.*XNT*TSTAR/TR-4.*XNT-4.*XNT*X/TR;
                    end
                end
                XNLD=(XNC-XNL)/TAU;
                YDD=YTDD-XNL;
                FLAG=1;
            end
            FLAG=0;
            Y=.5*(YOLD+Y+H*YD);
            YD=.5*(YDOLD+YD+H*YDD);
            XNL=.5*(XNLOLD+XNL+H*XNLD);
            S=S+H;
            if S>=(TS-.00001)
                S=0.;
                TGO=TF-T+.000001;
                RTM=VC*TGO;
                SIGNOISE=sqrt(SIGRIN^2+(SIGGL/RTM)^2+(SRN*RTM*RTM/(RA*RA))^2);
                SIGPOS=RTM*SIGNOISE;
                SIGN2=SIGPOS^2;
                RMAT=[SIGN2];
                M=PHI*P*PHI'+Q;
                K = M*HMAT'/(HMAT*M*HMAT' + RMAT);
                P = (IDNP - K*HMAT)*M;
                [YB,YDB,X1B,X2B,X3B,X4B]=PROJECT6S(T,TS,YH,YDH,XNL,...
                    X1H,X2H,X3H,X4H,W,H);
                XLAMNOISE=SIGNOISE*randn;
                YSTAR=RTM*(XLAM+XLAMNOISE);
                RES(1,1)=YSTAR-YB;
                YH=YB+K(1,1)*RES(1,1);
                YDH=YDB+K(2,1)*RES(1,1);
                X1H=X1B+K(3,1)*RES(1,1);
                X2H=X2B+K(4,1)*RES(1,1);
                X3H=X3B+K(5,1)*RES(1,1);
                X4H=X4B+K(6,1)*RES(1,1);
                XNTH=X1H+X3H;
                XNC=XNP*(YH+YDH*TGO+.5*XNTH*TGO*TGO)/(TGO*TGO);
                if XNC>XLIM
                    XNC=XLIM;
                end
                if XNC<-XLIM
                    XNC=-XLIM;
                end
                count=count+1;
                ArrayT(count)=T;
                ArrayYTDD(count)=YTDD;
                ArrayXNTH(count)=XNTH;
            end
        end
figure
plot(ArrayT',ArrayYTDD',ArrayT',ArrayXNTH'),grid
title('RMS miss for various flight times')
xlabel('Flight Time (S)')
ylabel('RMS MISS (Ft) ')
clc
output=[ArrayT',ArrayYTDD',ArrayXNTH'];
save datfil.txt output -ascii
disp('Simulation Complete')
