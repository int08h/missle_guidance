clear
count=0;
PZ=3.;
TR=.5;
% TR=.01;
XL=PZ/2.;
PI=3.14159;
W=2.*PI/PZ;
XNT=322.;
X=PZ/2.-TR;
ALF=PI*TR/(2.*XL);
T=0.;
for T=0.:.1:6.
	if T<PZ
		if T<TR/2.
            YTDD=2.*XNT*T/TR;
		elseif T<(TR/2.+X)
            YTDD=XNT;
		elseif T<(3.*TR/2.+X)
            YTDD=-2.*XNT*T/TR+2.*XNT+2.*XNT*X/TR;
		elseif T<(3.*TR/2.+2.*X)
            YTDD=-XNT;
		else
            YTDD=2.*XNT*T/TR-4.*XNT-4.*XNT*X/TR;
		end
	elseif T<2.*PZ
        TSTAR=T-PZ;
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
	elseif T<3.*PZ
        TSTAR=T-2.*PZ;
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
        TSTAR=T-3.*PZ;
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
    XNTH=4.*XNT*(sin(ALF)*sin(W*T)+sin(3.*ALF)*sin(3.*W*T)/9.)/(PI*ALF);
    count=count+1;
    ArrayT(count)=T;
    ArrayYTDD(count)=YTDD;
    ArrayXNTH(count)=XNTH;
end
figure
plot(ArrayT,ArrayYTDD,ArrayT,ArrayXNTH),grid
xlabel('Time (Sec)')
ylabel('Acceleration and Estimate (f/s^2)')
clc
output=[ArrayT',ArrayYTDD',ArrayXNTH'];
save datfil.txt output  -ascii
disp 'simulation finished'