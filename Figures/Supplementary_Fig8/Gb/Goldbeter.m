function dydt = Goldbeter(t,y)
% Y= [m p0 p1 p2 pN]

vs=0.76;
vm=0.65;
vd=0.95;
ks=0.38;
kt1=1.9;
kt2=1.3;
v1=3.2;
v2=1.58;
v3=5;
v4=2.5;
k1=1;
k2=1;
k3=2;
k4=2;
ki=1;
km1=0.5;
kd=0.2;
n=3;

dydt=zeros(5,1);
dydt(1) = vs / (1+ (y(5)/ki)^4) - vm * y(1) /(km1+ y(1));
dydt(2) = ks * y(1) - v1 * y(2) /(k1+ y(2))+ v2 * y(3) /(k2+ y(3));
dydt(3) = v1 * y(2) /(k1+ y(2)) - v2 * y(3) /(k2+ y(3)) -  v3 * y(3) /(k3+ y(3))  +v4 * y(4) /(k4+ y(4));
dydt(4) = v3 * y(3) /(k3+ y(3))  - v4 * y(4) /(k4+ y(4)) - kt1 * y(4) + kt2 * y(5) - vd * y(4) / (kd + y(4));
dydt(5) = kt1 * y(4) - kt2 * y(5);
