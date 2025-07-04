function s=earthdist(lon1,lat1,lon2,lat2)
% Compute the length in km of the geodesic line connecting two
% geographic locations on the earth.

rad=pi/180;
rearth=6371.001;

phi=lon1*rad;
lambda=lat1*rad;
x1=cos(lambda).*cos(phi);
y1=cos(lambda).*sin(phi);
z1=sin(lambda);

phi=lon2*rad;
lambda=lat2*rad;
x2=cos(lambda).*cos(phi);
y2=cos(lambda).*sin(phi);
z2=sin(lambda);

s=rearth*acos(x1.*x2+y1.*y2+z1.*z2);
