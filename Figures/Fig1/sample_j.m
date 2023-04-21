function dydt = sample_j(t,y,zt,z,yt,yy)

Z = interp1(zt,z,t);
Y = interp1(yt,yy,t);

%x y
dydt = zeros(1,1);
dydt(1) = Z;

end