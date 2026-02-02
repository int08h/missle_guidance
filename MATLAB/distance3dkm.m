function DISTKM=distance3d(XT,YT,ZT,XF,YF,ZF)
R=sqrt(XT^2+YT^2+ZT^2);
RF=sqrt(XF^2+YF^2+ZF^2);
A=2.0926E7;
CBETA=(XT*XF+YT*YF+ZT*ZF)/(R*RF);
if CBETA<=1.
	BETA=acos(CBETA);
	DISTKM=A*BETA/3280.;
else
	DISTKM=(XF-XT)/3280.;
end
