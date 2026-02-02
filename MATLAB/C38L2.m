clear
count=0;
RUN=100;
TSTART=0.;
P=3.;
TR=1.;
XL=P/2.;
VC=3000.;
XNT=161.;
TAU=.5;
XNP=3.;
TF=10.;
X=P/2.-TR;
for TF=.1:.1:10,
    Z1=0.;
    for JJ=1:RUN,
	    SUM=rand(1);
        TSTART=SUM*TF;
        Y=0.;
        YD=0.;
        XLAMH=0.;
        T=0.;
        H=.01;
        S=0.;
        while T <= (TF - 1e-5)
            YOLD=Y;
            YDOLD=YD;
            XLAMHOLD=XLAMH;
            STEP=1;
            FLAG=0;
            while STEP <=1
            if FLAG==1
                Y=Y+H*YD;
                YD=YD+H*YDD;
                XLAMH=XLAMH+H*XLAMHD;
                T=T+H;
                STEP=2;
            end;
            if T<TSTART
                YTDD=0.;
            elseif (T-TSTART)<P
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
            elseif (T-TSTART)<2.*P;
                TSTAR=(T-TSTART)-P;
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
            elseif (T-TSTART)<3.*P
                TSTAR=(T-TSTART)-2.*P;
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
                TSTAR=(T-TSTART)-3.*P;
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
            TGO=TF-T+.00001;
            XLAM=Y/(VC*TGO);
            XLAMHD=(XLAM-XLAMH)/TAU;
            XNC=XNP*VC*XLAMHD;
            YDD=YTDD-XNC;
            FLAG=1;
        end;
        FLAG=0;
        Y=.5*(YOLD+Y+H*YD);
        YD=.5*(YDOLD+YD+H*YDD);
        XLAMH=.5*(XLAMHOLD+XLAMH+H*XLAMHD);
        S=S+H;
    end
    Z(JJ)=Y;
    Z1=Z(JJ)+Z1;
    XMEAN=Z1/JJ;
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
