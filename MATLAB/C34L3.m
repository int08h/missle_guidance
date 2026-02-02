clear
count=0;
RUN=100;
ORDER=3;
IPOISSON=1;
%1=POISSON,2=NORMAL
VC=5.*3280.;
XNTIC=161.;
YIC=0.;
VM=3000.;
HEDEG=20.;
XNP=3.;
SIGNOISE=.0001;
TS=.01;
TAU=.2;
AMAXG=10.;
ICONSTANT=2;
% (1=POISSON,2=UNIF CONST,3=RANDOM SINE)
W=2.;
XNU=.5;
TSW=1;
AMAX=AMAXG*32.2;
BETA=XNTIC;
for TF=.1:.1:6.0,
    Z1=0;
    for JJ=1:RUN
        SUM=rand;
        TSTART=TF*SUM;
        SUM=rand;
        COEF=1.;
        PHASE=6.28*SUM;
        Y=0.;
        YD=0.;
        TS2=TS*TS;
        TS3=TS2*TS;
        TS4=TS3*TS;
        TS5=TS4*TS;
        PHI=zeros(ORDER);
        P=zeros(ORDER);
        QC=zeros(ORDER);
        IDNP=zeros(ORDER);
        F=zeros(ORDER);
        RTM=VC*TF;
        SIGPOS=RTM*SIGNOISE;
        SIGN2=SIGPOS^2;
        IDNP(1,1)=1.;
        IDNP(2,2)=1.;
        IDNP(3,3)=1.;
        HMAT(1,1)=1.;
        HMAT(1,2)=0.;
        HMAT(1,3)=0.;
        P(1,1)=SIGN2;
        P(2,2)=(VM*HEDEG/57.3)^2;
        P(3,3)=XNTIC*XNTIC;
        if IPOISSON==1
            F(1,2)=1.;
            F(2,3)=1.;
            F(3,3)=-2.*XNU;
            PHIN=4.*XNU*XNTIC*XNTIC;
        else IPOISSON==2
            F(1,2)=1.;;
            F(2,3)=1.;
            PHIN=XNTIC*XNTIC/10.;
        end
        QC(3,3)=PHIN;
        PHI=IDNP+F*TS+.5*F*F*TS*TS;
        Q=PHI*QC*PHI'*TS;
        T=0.;
        H=.001;
        S=0.;
        YH=0.;
        YDH=0.;
        XNTH=0.;
        XNC=0.;
        XNL=0.;
        BETA=XNTIC;
        QFIRST=1;
        SIG=1./sqrt(2.*XNU);
        XNOISE=randn;
        if XNOISE>0.
            XNTP=BETA;
        end
        if XNOISE<=0.
            XNTP=-BETA;
        end
        DELT=9999.;
        TNOW=0.;
        while T<=(TF-.0001)
            if QFIRST==1
                XNOISE1=SIG*randn;
                XNOISE2=SIG*randn;
                DELT=XNOISE1^2+XNOISE2^2;
                QFIRST=0.;
                TNOW=T;
            end
            if T>=(DELT+TNOW)
                XNTP=-XNTP;
                QFIRST=1;
            end
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
                if ICONSTANT==1
                    XNTC=XNTP;
                elseif ICONSTANT==2
                    if T<TSTART
                        XNTC=0.;
                    else
                        XNTC=COEF*XNTIC;
                    end
                else
                    if T<TSTART
                        XNTC=0.;
                    else
                        XNTC=XNTIC*sin(W*T+PHASE);
                    end
                end
                TGO=TF-T+.00001;
                RTM=VC*TGO;
                XLAM=Y/(VC*TGO);
                XNLD=(XNC-XNL)/TAU;
                YDD=XNTC-XNL;
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
                SIGPOS=RTM*SIGNOISE;
                SIGN2=SIGPOS^2;
                RMAT=[SIGN2];
                M=PHI*P*PHI'+Q;
                K = M*HMAT'/(HMAT*M*HMAT' + RMAT);
                P = (IDNP - K*HMAT)*M;
                XLAMNOISE=SIGNOISE*randn;
                YSTAR=RTM*(XLAM+XLAMNOISE);
                [YB,YDB,XNTB]=PROJECT34(T,TS,YH,YDH,XNL,...
                                   XNTH,XNU,IPOISSON);
                RES(1,1)=YSTAR-YB;
                YH=YB+K(1,1)*RES(1,1);
                YDH=YDB+K(2,1)*RES(1,1);
                XNTH=XNTB+K(3,1)*RES(1,1);
                X=TGO/TAU;
                TAUTGT=1./(2.*XNU);
                ZEM1H=YH+YDH*TGO-XNL*TAU*TAU*(exp(-X)+X-1.)+.5*XNTH*TGO*TGO;
                TOP=6.*X*X*(exp(-X)-1.+X);
                BOT1=2*X*X*X+3.+6.*X-6.*X*X;
                BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
                XNPP=TOP/(.0001+BOT1+BOT2);
                XNC=XNPP*ZEM1H/TGO^2;
                if XNC>AMAX
                    XNC=AMAX;
                end
                if XNC<-AMAX
                    XNC=-AMAX;
                end
            end
        end
        Z(JJ)=Y;
        Z1=Z(JJ)+Z1;
        XMEAN=Z1/JJ;
    end;
    SIGMA=0;
    Z1=0;
    Z2=0.;
    for JJ=1:RUN,
        Z1=(Z(JJ)-XMEAN)^2+Z1;
        Z2=Z(JJ)^2+Z2;
        if JJ==1,
            SIGMA=0;
            RMS=0.;
        else
            SIGMA=sqrt(Z1/(JJ-1));
            RMS=sqrt(Z2/(JJ-1));
        end;
    end;
    count=count+1;
    ArrayTF(count)=TF;
    ArraySIGMA(count)=SIGMA;
    ArrayXMEAN(count)=XMEAN;
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



        