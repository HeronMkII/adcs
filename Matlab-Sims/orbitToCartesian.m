function [P] = orbitToCartesian(K)
%input units must be in meters and degrees
%output units are in meters and seconds

u=K(1);
a=K(2);
e=K(3);
i=K(4);
w=K(5);
ohm=K(6);
M=K(7);

i=i/180*pi;
w=w/180*pi;
ohm=ohm/180*pi;
M=M/180*pi;

E=M;
E_old=0;
while abs((E-E_old))>10e-4
    E_old=E;
    E=E_old-(E_old-e*sin(E_old)-M)/(1-e*cos(E_old));
end

v=2*atan2((sqrt(1+e)*sin(E/2)),(sqrt(1-e)*cos(E/2)));

r=a*(1-e*cos(E));

o=r*[cos(v);sin(v);0];
odot=(sqrt(u*a)/r)*[-sin(E);sqrt(1-e^2)*cos(E);0];

P(1:3,1)=[o(1)*(cos(w)*cos(ohm)-sin(w)*cos(i)*sin(ohm))-o(2)*(sin(w)*cos(ohm)+cos(w)*cos(i)*sin(ohm));
          o(1)*(cos(w)*sin(ohm)+sin(w)*cos(i)*cos(ohm))+o(2)*(cos(w)*cos(i)*cos(ohm)-sin(w)*sin(ohm));
          o(1)*sin(w)*sin(i)+o(2)*cos(w)*sin(i)];
P(4:6,1)=[odot(1)*(cos(w)*cos(ohm)-sin(w)*cos(i)*sin(ohm))-odot(2)*(sin(w)*cos(ohm)+cos(w)*cos(i)*sin(ohm));
          odot(1)*(cos(w)*sin(ohm)+sin(w)*cos(i)*cos(ohm))+odot(2)*(cos(w)*cos(i)*cos(ohm)-sin(w)*sin(ohm));
          odot(1)*sin(w)*sin(i)+odot(2)*cos(w)*sin(i)];
      
end

