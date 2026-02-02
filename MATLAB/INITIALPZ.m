function [rt1f,rt2f,tfdes]=initialpz(rt2des,rt1ic,rt2ic,...
					vt1ic,vt2ic,beta)
rt1=rt1ic;
rt2=rt2ic;
vt1=vt1ic;
vt2=vt2ic;
t=0.;
h=.01;
while rt2>rt2des
	rt1old=rt1;
	rt2old=rt2;
	vt1old=vt1;
	vt2old=vt2;
 	step=1;
 	flag=0;
	while step <=1
		if flag==1
			rt1=rt1+h*vt1;
			rt2=rt2+h*vt2;
			vt1=vt1+h*at1;
			vt2=vt2+h*at2;
 			t=t+h;
 			step=2;
 		end
 		if rt2<=30000.
			rho=.002378*exp(-rt2/30000.);
		else
			rho=.0034*exp(-rt2/22000.);
		end
		vt=sqrt(vt1^2+vt2^2);
		q=.5*rho*vt^2;
		gamt=atan2(-vt2,vt1);
		at1=-32.2*q*cos(gamt)/beta;
		at2=-32.2+32.2*q*sin(gamt)/beta;
		flag=1;
	end;
	flag=0;
	rt1=.5*(rt1old+rt1+h*vt1);
	rt2=.5*(rt2old+rt2+h*vt2);
	vt1=.5*(vt1old+vt1+h*at1);
	vt2=.5*(vt2old+vt2+h*at2);
end
rt1f=rt1;
rt2f=rt2;
tfdes=t;
	
