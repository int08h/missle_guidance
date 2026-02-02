function [xtf,ytf]=predict45(tp,xtp,ytp,xtdp,ytdp,tf,tftot,tupt,xf,yf,itgt)
if itgt==1
	tpz=180;
else
	tpz=240;
end
t=tp;
switch1=0;
xt=xtp;
yt=ytp;
zt=0.;
xtd=xtdp;
ytd=ytdp;
ztd=0.;
zf=0.;
a=2.0926E7;
gm=1.4077E16;
qboost=1;
h=.01;
s=0.;
axt=0.;
ayt=0.;
ztd=0;
while t<=(tf-.00001)
	xtold=xt;
	ytold=yt;
	xtdold=xtd;
	ytdold=ytd;
	step=1;
	flag=0;
	while step <=1
		if flag==1
			xt=xt+h*xtd;
			yt=yt+h*ytd;
			xtd=xtd+h*xtdd;
			ytd=ytd+h*ytdd;
			t=t+h;
			step=2;
		end
		tembot=(xt^2+yt^2)^1.5;
		xtdd=-gm*xt/tembot+axt;
		ytdd=-gm*yt/tembot+ayt;
		if itgt==1
			if t<180.
				wgt=-212.*t+44000.;
				trst=54100.;
			else
				wgt=3300.;
				trst=0.;
			end
		else
			if t<120
				wgt=-2622*t+440660.;
				trst=725850.;
			elseif t<240.
				wgt=-642.*t+168120.;
				trst=182250.;
			else
				wgt=5500.;
				trst=0.;
			end
		end
		atp=32.2*trst/wgt;
		flag=1;
	end;
	flag=0;
 	xt=(xtold+xt)/2+.5*h*xtd;
	yt=(ytold+yt)/2+.5*h*ytd;
	xtd=(xtdold+xtd)/2+.5*h*xtdd;
	ytd=(ytdold+ytd)/2+.5*h*ytdd;
	if qboost==1
		tgolam=tftot-t;
		[vrx,vry,vrz]=LAMBERT3D(xt,yt,zt,tgolam,xf,yf,zf,switch1);
		vtx=vrx;
		vty=vry;
		delvxt=vtx-xtd;
		delvyt=vty-ytd;
		delvelt=sqrt(delvxt^2+delvyt^2);
		if (t<tpz & delvelt>500.)
			axt=atp*delvxt/delvelt;
			ayt=atp*delvyt/delvelt;
		elseif delvelt<500.
			trst=0.;
			qboost=1;
			axt=0.;
			ayt=0.;
			xtd=vtx;
			xtdold=xtd;
			ytd=vty;
			ytdold=ytd;
		else
			qboost=0;
			qoomph=0;
			axt=0.;
			ayt=0.;
		end
		if t<tupt
			rtmag=sqrt(xt^2+yt^2);
			axt=atp*xt/rtmag;
			ayt=atp*yt/rtmag;
		end
	end
end
xtf=xt;
ytf=yt;