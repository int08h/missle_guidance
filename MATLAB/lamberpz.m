function [vrx,vry]=lambert(xic,yic,tfdes,xf,yf,xlongm,xlongt)
a=2.0926e7;
gm=1.4077e16;
ric=sqrt(xic^2+yic^2);
rf=sqrt(xf^2+yf^2);
cphi=(xic*xf+yic*yf)/(ric*rf);
phi=acos(cphi);
sphi=sin(phi);
r0=ric;
degrad=360./(2.*pi);
icount=0;
gmin=atan2((sphi-sqrt(2.*r0*(1.-cphi)/rf)),(1-cphi));
gmax=atan2((sphi+sqrt(2.*r0*(1.-cphi)/rf)),(1-cphi));
gam=(gmin+gmax)/2.;
tf=0;
while ~(abs(tfdes-tf)<=(.00000001*tfdes))
	top=gm*(1.-cos(phi));
	temp=r0*cos(gam)/rf-cos(phi+gam);
	bot=r0*cos(gam)*temp;
	v=sqrt(top/bot);
	if  xlongt>xlongm
	  	vrx=v*cos(pi/2.-gam+xlongm);
	  	vry=v*sin(pi/2.-gam+xlongm);
	else
	  	vrx=v*cos(-pi/2.+gam+xlongm);
	  	vry=v*sin(-pi/2.+gam+xlongm);
	end
	xlam=r0*v*v/gm;
	top1=tan(gam)*(1-cos(phi))+(1-xlam)*sin(phi);
	bot1p=(1-cos(phi))/(xlam*cos(gam)*cos(gam));
	bot1=(2-xlam)*(bot1p+cos(gam+phi)/cos(gam));
	top2=2*cos(gam);
	bot2=xlam*((2/xlam-1)^1.5);
	top3=sqrt(2/xlam-1);
	bot3=cos(gam)/tan(phi/2)-sin(gam);
	temp=(top2/bot2)*atan2(top3,bot3);
	tf=r0*(top1/bot1+temp)/(v*cos(gam));
	icount=icount+1;
	if tf>tfdes
		gmax=gam;
	else
		gmin=gam;
	end
	if icount==1
		xnext=(gmax+gmin)/2.;
	else
		xnext=gam+(gam-gold)*(tfdes-tf)/(tf-told);
		if (xnext>gmax|xnext<gmin)
			xnext=(gmax+gmin)/2.;
		end
	end
	gold=gam;
	told=tf;
	gam=xnext;
end
