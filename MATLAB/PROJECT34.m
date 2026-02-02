function [YB,YDB,XNTB]=PROJECT34(TP,TS,YHP,YDHP,XNLP,XNTHP,XNU,...
    IPOISSON)
T=0.;
YH=YHP;
YDH=YDHP;
XNT=XNTHP;
XNL=XNLP;
H=.001;
while T<=(TS-.0001)
    if IPOISSON==1
        XNTD=-2.*XNU*XNT;
        XNT=XNT+H*XNTD;
        YDDH=XNT-XNL;
        YDH=YDH+H*YDDH;
        YH=YH+H*YDH;
    else
        YDDH=XNT-XNL;
        YDH=YDH+H*YDDH;
        YH=YH+H*YDH;
    end
    T=T+H;
end
YB=YH;
YDB=YDH;
XNTB=XNT;
