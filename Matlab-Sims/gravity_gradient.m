function [Tg] = gravity_gradient(u, r, Rbi, I)
%Computes the gravity gradient torque on a body.
%u is the gravitational parameter
%r is the position vector to the body in ECI grame
%Rbi is the rotation matrix from ECI to Body Frame
%I is the moment of inertia tensor of the body

nadir = Rbi*(-r/norm(r));   %Nadir in Body frame

Tg = (3*u/(norm(r)^3))*skew(nadir)*I*nadir;

end

