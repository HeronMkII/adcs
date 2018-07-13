%% ATTITUDE AND ORBIT PROPAGATION %%

%Using a dipole model for the Earth's magnetic field, this program
%propagates the attitude of a passive magnetic attitude controlled
%satellite in orbit. The main restoring torque comes from a permanent
%magnet along the z-axis (minor axis), while the damping torque comes from
%hysteresis rods along the x and y axes (major axes). The gravity gradient
%torque is included as a disturbance torque.

clear all;
makeAttitudeVideo = false;

%Universal Constants
u=3.986004418e14;               %Gravitational parameter of the Earth
mu=4e-7*pi;                     %Relative permeability of free space

%Parameter Constants
days=7;                        %Simulation time in days
initialDay=334;                   %Initial day of the year (day 0 is at vernal equinox)
rodsPerAxis=1.4;                 %Number of hysteresis rods per axis (x and y)
% Hc=12;                          %Constant coercivity value of hyst rods
% Bs=0.025;                       %Magnetization at saturation of hyst rods
% Br=0.004;                       %Magnetic remanence of hysteresis rods
% %https://digitalcommons.usu.edu/cgi/viewcontent.cgi?referer= ...
% %%https://www.google.ca/&httpsredir=1&article=1230&context=smallsat

%HyMu80
Hc=1.59;                           %Constant coercivity value of hyst rods
Bs=0.73;                        %Magnetization at saturation of hyst rods
Br=0.35;                       %Magnetic remanence of hysteresis rods

M=[0;0;2];                    %Magnetic moment of permanent magnet

mass=3.138; %kg                       %Satellite mass and dimentions
width=0.1; %m
height=0.3; %m

rodLength=0.07; %m               %Hysteresis rod dimensions
rodRadius=0.0005; %m

%hot case - 31 JAN

% majorAxis=6928.14e3;      %m        %Orbital parameters in keplerian elements
% eccentricity=0;        
% inclination=97.5897; %deg
% longOfAscendingNode=98.9563; %deg
% ArgumentOfPeriapsis=0;
% MeanAnomalyEpoch=359.891;

%cold case - NOV 1

majorAxis=6928.14e3;      %m        %Orbital parameters in keplerian elements
eccentricity=0;        
inclination=97.656; %deg
longOfAscendingNode=16.4596; %deg
ArgumentOfPeriapsis=0;
MeanAnomalyEpoch=359.908;

% er = 0.362767142;           %for 10:30 SSO @ 500km, Nov 2019 average
% er = 0.3;                      %hot case, jan 31 SSO 550 km
er = 0.371;

period = 5746; %in seconds, for SSO @ 550 km

earthInclination=23.5;        %Inclination of Earth orbit plane in degrees

%Calculated constants
time=days*3600*24;                          %Simulation time in seconds
solarPanelArea=6*26.62e-4;                  %26.62cm2/panel, 6 per face
% I=[(1/12)*mass*(width^2+height^2),0,0;...   %Inertia tensor assuming
%     0,(1/12)*mass*(width^2+height^2),0;...     %uniform mass distribution
%     0,0,(1/12)*mass*(width^2+width^2)];
I=[0.026066,0,0;...   %Inertia tensor assuming
    0,0.026066,0;...     %uniform mass distribution
    0,0,0.00523];
rodVol=pi*rodRadius^2*rodLength;            %Hysteresis rod volume
Hr=cot((pi*Br)/(2*Bs))*Hc;                  %Remanence constant

%Initial Conditions
kepElem=[u;majorAxis;...            %Orbit parameters in keplerian elements
         eccentricity;...
         inclination;...
         longOfAscendingNode;...
         ArgumentOfPeriapsis;...
         MeanAnomalyEpoch];
p0=orbitToCartesian(kepElem);       %Convert keplerian elements to cartesian position and velocity
r0=p0(1:3,1);                       %Initial position in ECI (or ECF because initially aligned)
v0=p0(4:6,1);                       %Initial velocity in ECI (or ECF because initially aligned)

Ric = [1 0 0;...                            %Rotation cosine matrix from  
       0 cosd(earthInclination) -sind(earthInclination);... %celestial inertial to
       0 sind(earthInclination) cosd(earthInclination)];    %equatorial inertial

w0=[0.2;0.2;0.1];                %Initial angular velocity
% w0=[0.1;0.1;0.01];
q0=[1;0;0;0];                       %Initial quaternion relating satellite attitude to ECI frame

%Determine initial magnetization of rods assuming they started out at B=0
%and have reached equilibrium (they have been in the constant external
%magnetic field for a long time)
Rbi0 = ECItoBody(q0);
Hb0 = Rbi0*(dipole_magstrength(r0));   %ECF and ECI are initially aligned, so no need to rotate
HxRod0=Hb0(1);
HyRod0=Hb0(2);
if HxRod0 < 0
    x = atan(-(Hc-HxRod0)/Hr);
    BxRod0 = x*2*Bs/pi;
else
    x = atan((HxRod0+Hc)/Hr);
    BxRod0 = x*2*Bs/pi;
end
if HyRod0 < 0
    y = atan(-(Hc-HyRod0)/Hr);
    ByRod0 = y*2*Bs/pi;
else
    y = atan((HyRod0+Hc)/Hr);
    ByRod0 = y*2*Bs/pi;
end

%setup eclipse array
orbits = floor(time/period); %# of orbits
edur = floor(period*er); %eclipse duration
ec = zeros(1,time+1); %1 for eclipse 0 for not
for k=1:orbits
    for c= 1:edur
        counter = k*period-edur + c;
        ec(counter)=1;        
    end
end

%INTEGRATE
options = odeset('RelTol',1e-6,'AbsTol',1e-7,'OutputFcn',@odetpbar);
[t,A] = ode23(@(t,x)propagate(t, x, mu, u, M, rodVol, I, Hc, Hr, Bs, rodsPerAxis,Ric,ec,time,orbits,edur,period),...
              0:time, [w0;q0;r0;v0;BxRod0;ByRod0], options);

%Integration outputs
w = A(:,1:3);
q = A(:,4:7);
r = A(:,8:10);
v = A(:,11:13);
BxRod = A(:,14);
ByRod = A(:,15);

%%GENERATE RESULTS

rEarth = zeros(3,time+1);   %%Position in ECF
vEarth = zeros(3,time+1);   %%Velocity in ECF

He = zeros(3,time+1);       %%Magnetic field strength in ECF
Hb = zeros(3,time+1);       %%Magnetic field strength in Body Frame
Be = zeros(3,time+1);       %%Magnetic field flux density in ECF
Bb = zeros(3,time+1);       %%Magnetic field flux density in Body Frame
HxRod = zeros(1,time+1);    %%Magnetic field strength along an x-axis rod
HyRod = zeros(1,time+1);    %%Magnetic field strength along a y-axis rod

Tmag = zeros(3,time+1);     %%Torque due to permanent magnet
TxRod = zeros(3,time+1);    %%Torque due to all x-axis rods
TyRod = zeros(3,time+1);    %%Torque due to all y-axis rods
Tg = zeros(3,time+1);       %%Torque due to the gravity gradient

Error = zeros(1,time+1);    %%Angle between minor (z) axis and external magnetic field

timeDays = zeros(time+1,1); %%Time of simulation in days

xNegFlux = zeros(1,time+1); %%Flux on solar panel surfaces
xPosFlux = zeros(1,time+1);
yNegFlux = zeros(1,time+1);
yPosFlux = zeros(1,time+1);

totalFlux = zeros(1,time+1);
%%For-loop to fill matrices at every time step
for k=1:time+1      %%Each index of k corresponds to the current time+1 (k=1->t=0)
    
    %%Declare rotation cosine matrices
    Rei = ECItoECF(k-1);
    ReiDot = dECItoECFdt(k-1);
    Rbi = ECItoBody(q(k,:));
    
    %%Fill rEarth and vEarth
    rEarth(:,k) = Rei*r(k,:)';
    vEarth(:,k) = ReiDot*r(k,:)' + Rei*v(k,:)';
    
    %%Fill magnetic field elements
    He(:,k) = dipole_magstrength(rEarth(:,k));
    Be(:,k) = dipole_magfield(rEarth(:,k));
    Hb(:,k) = Rbi*Rei'*He(:,k);
    Bb(:,k) = Rbi*Rei'*Be(:,k);
    HxRod(k) = Hb(1,k);
    HyRod(k) = Hb(2,k);
    
    %%Fill torque elements
    Tmag(:,k) = skew(M)*(Bb(:,k));
    TxRod(:,k) = rodsPerAxis * skew([BxRod(k)*rodVol;0;0])*Hb(:,k);
    TyRod(:,k) = rodsPerAxis * skew([0;ByRod(k)*rodVol;0])*Hb(:,k);
    Tg(:,k) = gravity_gradient(u,r(k,:)',Rbi,I);
    
    %Fill Error matrix
    Error(k) = angle(M,Bb(:,k));
    
    %Fill TimeDays matrix
    timeDays(k) = t(k)/3600/24;
    
    %%Find the solar vector
    theta = ((((k-1)/3600/24)+initialDay)/365.25)*2*pi;  %Angle Earth has orbited Sun in rad
    rEtoS_C = [cos(theta);sin(theta);0];    %Unit vector from Earth to Sun in celestial inertial 
    rEtoS_I = Ric*rEtoS_C;                  %In equitorial inertial
    rEtoS_B = Rbi*rEtoS_I;                  %In body frame
    rEtoS_Bn = rEtoS_B / sqrt(rEtoS_B' * rEtoS_B); %normalize, make unit vector
    
    %%Fill Flux matrices
    xNegFlux(k) = 1385*[-1 0 0]*rEtoS_Bn; %1395 - avg solar irridiance of Nov 2019 10:30
    xPosFlux(k) = 1385*[1 0 0]*rEtoS_Bn;
    yNegFlux(k) = 1385*[0 -1 0]*rEtoS_Bn;
    yPosFlux(k) = 1385*[0 1 0]*rEtoS_Bn;
    
    %%Disregard Flux if it is negative (surface is not facing the sun)
    if xNegFlux(k)<0
        xNegFlux(k)=0;
    end
    if xPosFlux(k)<0
        xPosFlux(k)=0;
    end
    if yNegFlux(k)<0
        yNegFlux(k)=0;
    end
    if yPosFlux(k)<0
        yPosFlux(k)=0;
    end    
end

%%%%%%%%ECLIPSE
orbits = floor(time/period);
edur = floor(period*er); %eclipse duration
for k=1:orbits
    for c= 1:edur
        counter = k*period-edur + c;
        xNegFlux(counter)=0;
        xPosFlux(counter)=0;
        yNegFlux(counter)=0;
        yPosFlux(counter)=0;
    end
end

%%GRAPH RESULTS

%%Easily designate which results to graph
angVel = true;
errorAng = true;
torque = false;
hystCurves = false;
inertialOrbit = false;
earthOrbit = false;
qNorm = false;
tempg = false;
flux = false;

if angVel==true
    figure;
    plot(timeDays,w(:,1),'r'); hold on;
    plot(timeDays,w(:,2),'b')
    plot(timeDays,w(:,3),'k')
    title('Body-referenced angular velocities','FontSize',14);
    xlabel('time (days)','FontSize',12);
    ylabel('rad/s)','FontSize',12);
    legend('w(x)','w(y)','w(z)');
end

if errorAng==true
    figure;
    plot(timeDays,Error)
    title('Error angle','FontSize',14);
    xlabel('time (days)','FontSize',12);
    ylabel('degrees','FontSize',12);
end

if torque==true
    figure;
    plot(timeDays,sqrt(Tmag(1,:).^2+Tmag(2,:).^2+Tmag(3,:).^2),'r'); hold on;
    plot(timeDays,sqrt(TxRod(1,:).^2+TxRod(2,:).^2+TxRod(3,:).^2),'b');
    plot(timeDays,sqrt(TyRod(1,:).^2+TyRod(2,:).^2+TyRod(3,:).^2),'k');
    plot(timeDays,sqrt(Tg(1,:).^2+Tg(2,:).^2+Tg(3,:).^2),'g');
    title('Magnitudes of Torques on Satellite','FontSize',14);
    xlabel('time (days)','FontSize',12);
    ylabel('N*m','FontSize',12);
    legend('Permanent Magnet','x-Axis Rods','y-Axis Rods','Gravity Gradient');
end

if hystCurves==true
    figure;
    plot(HxRod,BxRod);
    title('Hysteresis loop of x-axis rod','FontSize',14);
    xlabel('H','FontSize',12);
    ylabel('B','FontSize',12);
    figure;
    plot(HyRod,ByRod);
    title('Hysteresis loop of y-axis rod','FontSize',14);
    xlabel('H','FontSize',12);
    ylabel('B','FontSize',12);
end

if inertialOrbit==true
    figure;
    plot3(r(:,1),r(:,2),r(:,3),'k');
    title('Orbit in Inertial Frame (m)','FontSize',12);
end

if earthOrbit==true
    figure;
    plot3(rEarth(1,:),rEarth(2,:),rEarth(3,:),'k');
    title('Orbit in Earth-Fixed Frame (m)','FontSize',12);
end

if qNorm==true
    figure;
    plot(timeDays,sqrt(q(:,1).^2+q(:,2).^2+q(:,3).^2+q(:,4).^2)) %normalization over time
end

if flux==true
    figure; 
    hold on;
    plot(t,xPosFlux)
    plot(t,xNegFlux)
    plot(t,yPosFlux)
    plot(t,yNegFlux)
    legend('xPos','xNeg','yPos','yNeg');
    title('Flux On Solar Panels','FontSize',14);
    xlabel('time (s)','FontSize',12);
    ylabel('Flux (W/m^2)','FontSize',12);
end


%Make the attitude video
if makeAttitudeVideo==true
    attitudeVideo = makeAnimation(q,Bb,time);
end