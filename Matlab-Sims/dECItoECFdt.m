function [ReiDot] = dECItoECFdt(t)
%Given a certain time t since epoch (at which ECF=ECI),
%returns the time derivative of the 
%rotation cosine matrix from ECI to ECF.

wEarth=7.292115090e-5;

ReiDot = [-wEarth*sin(wEarth*t), -wEarth*cos(wEarth*t), 0;...
          wEarth*cos(wEarth*t), -wEarth*sin(wEarth*t), 0;...
          0, 0, 0];

end

