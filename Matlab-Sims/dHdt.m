function [HeDot] = dHdt(r,v)
%Returns magnetic field strength time derivative in 
%ECF coordinates, assuming a dipole model for the
%Earth. Core is considered a giant magnet.

Me=[0;0;-7.94e22];

HeDot = (3/(4*pi))*(((v*Me'*r+r*Me'*v)*norm(r)^2 - ...
         5*r*Me'*r*r'*v)/norm(r)^7 ...
         - Me*r'*v/norm(r)^5);

end

