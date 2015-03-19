function [a] = dengue_no_kids_events ()

clear;

t(1)=0.;
k=1;
N=100.;
i1(1)=5.; #first infection, serotype1
i2(1)=3.; #first infection, serotype2
r1(1)=0.; #immune to serotype1
r2(1)=0.; #immune to serotype2
i21(1)=0; #second infection, serotype2
i12(1)=0; #second infection, serotype1
r(1)=0.; #immune to both of them

s(1)=N-i1(1)-i2(1)-r1(1)-r2(1)-i21(1)-i12(1)-r(1);


a(1,1:9)=0;
a(1,1)=t(1);

events=400;
Ro=3;
gam1=.7;
gam2=.7;
mu=.00024;
nu=.015;
bet1=Ro*(mu+gam1); 
bet2=10*bet1;





while i1(k)+i2(k)+i12(k)+i21(k)>0.& k<events;
	u1=rand; u2=rand;
	ev1=(bet1*i1(k)+bet2*i12(k))*s(k)/N;
	ev2=(bet1*i2(k)+bet2*i21(k))*s(k)/N;
	ev3=gam1*i1(k);
	ev4=gam1*i2(k);
	ev5=r1(k)*(bet1*i2(k)+bet2*i21(k))/N;
	ev6=r2(k)*(bet1*i1(k)+bet2*i12(k))/N;
	ev7=gam1*i21(k);
	ev8=gam1*i12(k);
	ev9=mu*i1(k);
	ev10=mu*i2(k);
	ev11=mu*r1(k);
	ev12=mu*r2(k);
	ev13=(mu+nu)*i21(k);
	ev14=(mu+nu)*i12(k);
	ev15=mu*r(k) +ev14;
	tev2=ev2+ev1;
	tev3=tev2+ev3;
	tev4=tev3+ev4;
	tev5=tev4+ev5;
	tev6=tev5+ev6;
	tev7=tev6+ev7;
	tev8=tev7+ev8;
	tev9=tev8+ev9;
	tev10=tev9+ev10;
	tev11=tev10+ev11;
	tev12=tev11+ev12;
	tev13=tev12+ev13;
	tev14=tev13+ev14;
	tev15=tev14+ev15;
	t(k+1)=max(-log(u1)/tev15,0)+t(k);
	if u2*tev15<=ev1 & ev1!=0
		s(k)=s(k)-1;
		i1(k)=i1(k)+1;
		a(k,2)=-1;
		a(k,3)=1;
	elseif ev1< u2*tev15 & u2*tev15 <= tev2 & ev2 !=0
		s(k)=s(k)-1;
		i2(k)=i2(k)+1;
		a(k,2)=-1;
		a(k,4)=1;
	elseif tev2< u2*tev15 & u2*tev15 <= tev3 & ev3 !=0
		i1(k)=i1(k)-1;
		r1(k)=r1(k)+1;
		a(k,3)=-1;
		a(k,5)=1;
	elseif tev3< u2*tev15 & u2*tev15 <= tev4 & ev4 !=0
		i2(k)=i2(k)-1;
		r2(k)=r2(k)+1;
		a(k,4)=-1;
		a(k,6)=1;
	elseif tev4< u2*tev15 & u2*tev15 <= tev5 & ev5 !=0
		r1(k)=r1(k)-1;
		i21(k)=i21(k)+1;
		a(k,5)=-1;
		a(k,7)=1;
	elseif tev5< u2*tev15 & u2*tev15 <= tev6  & ev6 !=0
		r2(k)=r2(k)-1;
		i12(k)=i12(k)+1;
		a(k,6)=-1;
		a(k,8)=1;
	elseif tev6< u2*tev15 & u2*tev15 <= tev7 & ev7 !=0
		i21(k)=i21(k)-1;
		r(k)=r(k)+1;
		a(k,7)=-1;
		a(k,9)=1;
	elseif tev7< u2*tev15 & u2*tev15 <= tev8 & ev8 !=0
		i12(k)=i12(k)-1;
		r(k)=r(k)+1;
		a(k,8)=-1;
		a(k,9)=1;
	elseif tev8< u2*tev15 & u2*tev15 <= tev9 & ev9 !=0
		i1(k)=i1(k)-1;
		s(k)=s(k)+1;
		a(k,3)=-1;
		a(k,2)=1;
	elseif tev9< u2*tev15 & u2*tev15 <= tev10 & ev10 !=0
		i2(k)=i2(k)-1;
		s(k)=s(k)+1;
		a(k,4)=-1;
		a(k,2)=1;
	elseif tev10< u2*tev15 & u2*tev15 <= tev11 & ev11 !=0
		r1(k)=r1(k)-1;
		s(k)=s(k)+1;
		a(k,5)=-1;
		a(k,2)=1;
	elseif tev11< u2*tev15 & u2*tev15 <= tev12 & ev12 !=0
		r2(k)=r2(k)-1;
		s(k)=s(k)+1;
		a(k,6)=-1;
		a(k,2)=1;
	elseif tev12< u2*tev15 & u2*tev15 <= tev13 & ev13 !=0
		i21(k)=i21(k)-1;
		s(k)=s(k)+1;
		a(k,7)=-1;
		a(k,2)=1;
	elseif tev13< u2*tev15 & u2*tev15 <= tev14 & ev14 !=0
		i12(k)=i12(k)-1;
		s(k)=s(k)+1;
		a(k,8)=-1;
		a(k,2)=1;
	elseif tev14< u2*tev15 & u2*tev15 <= tev15 & ev15 !=0
		r(k)=r(k)-1;
		s(k)=s(k)+1;
		a(k,9)=-1;
		a(k,2)=1;
	end
	s(k+1)=s(k);
	i1(k+1)=i1(k);
	i2(k+1)=i2(k);
	r1(k+1)=r1(k);
	r2(k+1)=r2(k);
	i21(k+1)=i21(k);
	i12(k+1)=i12(k);
	r(k+1)=r(k);
	a(k,1)=t(k);
	k=k+1;


end




endfunction


