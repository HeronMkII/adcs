function [Rei] = ECItoECF(t)
%Given a certain time t since epoch (at which ECF=ECI),
%returns the rotation cosine matrix from ECI to ECF.

wEarth=7.292115090e-5;

Rei = [cos(wEarth*t) sin(wEarth*t) 0;...
       -sin(wEarth*t) cos(wEarth*t) 0;...
       0 0 1];

end

