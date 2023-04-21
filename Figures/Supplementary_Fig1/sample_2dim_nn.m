function dydt = sample_2dim_nn(t,y,zt,z,yt,yy)

Z = interp1(zt,z,t);
Y = interp1(yt,yy,t);

%x y
dydt = zeros(1,1);
dydt(1) = -Z-Y;

end