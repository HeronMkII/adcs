function [Y] = propagate(t, x, mu, u, M, rodVol, I, Hc, Hr, Bs, rodsPerAxis)
%ODE time propagation of the:
%attitude
%angular velocity
%orbital position
%orbital velocity
%rod magnetizations

%Seperate state matrix x into its components
w = x(1:3,1);
q = x(4:7,1);
r = x(8:10,1);
v = x(11:13,1);
BxRod = x(14,1);
ByRod = x(15,1);

%Cosine rotation matrices
Rbi = ECItoBody(q);
Rei = ECItoECF(t);
ReiDot = dECItoECFdt(t);

%Orbital position and velocity in ECF coordinates
rEarth = Rei*r;
vEarth = ReiDot*r + Rei*v;

%Magnetic field strength along x and y rods
He = dipole_magstrength(rEarth);
Hb = Rbi*Rei'*He;
HxRod = Hb(1);
HyRod = Hb(2);

%Magnetic field strength derivative in body frame
HeDot = dHdt(rEarth,vEarth);
HbDot = (-skew(w)*Rbi*Rei' + Rbi*ReiDot')*He + Rbi*Rei'*HeDot;
HDotxRod = HbDot(1);
HDotyRod = HbDot(2);

%Permanent magnet torque
%T = M X B
%B = mu*H
%--> T = M X (mu*H)
Tmag = skew(M)*(mu*Hb);

%Rod torques
%T = M X B
%M = (Brod*rodVolume)/mu
%B = mu*H
%--> T = (Brod*rodVolume) X Hb
TxRod = skew([BxRod*rodVol;0;0])*Hb;
TyRod = skew([0;ByRod*rodVol;0])*Hb;

%Gravity gradient torque
Tg = gravity_gradient(u,r,Rbi,I);

%Total torque
Ttotal = Tmag + rodsPerAxis*(TxRod + TyRod) + 25.2*Tg;

%Integration
Y(1:3,1) = I\(Ttotal-skew(w)*(I*w));                        %Gives w (integral of angular acceleration)
Y(4:7,1) = 0.5*S(w)*q;                                      %Gives q (qDot = 0.5*S(w)*q)
Y(8:10,1) = v;                                              %Gives r
Y(11:13,1) = -u*r/norm(r)^3;                                %Gives v (integral of acceleration)
if HDotxRod < 0                                             %Integration formula for magnetization
    Y(14,1) = (2*Bs)/(Hr*pi)*...                                %of hysteresis rods depends on sine of
              (((Hc-HxRod)*cos((pi*BxRod)/(2*Bs))+...           %the derivative of H along the rod
               Hr*sin((pi*BxRod)/(2*Bs))) / (2*Hc))^2*...
              HDotxRod;
else
    Y(14,1) = (2*Bs)/(Hr*pi)*...
              (((Hc+HxRod)*cos((pi*BxRod)/(2*Bs))-...
               Hr*sin((pi*BxRod)/(2*Bs))) / (2*Hc))^2*...
              HDotxRod;
end
if HDotyRod < 0
    Y(15,1) = (2*Bs)/(Hr*pi)*...
              (((Hc-HyRod)*cos((pi*ByRod)/(2*Bs))+...
               Hr*sin((pi*ByRod)/(2*Bs))) / (2*Hc))^2*...
              HDotyRod;
else
    Y(15,1) = (2*Bs)/(Hr*pi)*...
              (((Hc+HyRod)*cos((pi*ByRod)/(2*Bs))-...
               Hr*sin((pi*ByRod)/(2*Bs))) / (2*Hc))^2*...
              HDotyRod;
end

end