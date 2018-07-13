timee = length(xPosFlux) -1;
qsx = xPosFlux; %W/m^2
qsxx = xNegFlux;
qsy = yPosFlux;
qsyy = yNegFlux;
alpha = 0.92;
epsilon = 0.85;
At = solarPanelArea; %m^2
Asx = zeros(1,timee+1); %m^2
% Asxx = zeros(1,timee+1);
% Asy = zeros(1,timee+1);
% Asyy = zeros(1,timee+1);

for i=1:timee+1
    Asx(i) = At*qsx(i)/1395;
%     Asxx(i) = At*qsxx(i)/1395;
%     Asy(i) = At*qsy(i)/1395;
%     Asyy(i) = At*qsyy(i)/1395;
end

fa = 0.227; %????
beta = 20.336*pi/180; %rad
qe = 218; %?????
bolt = 5.67e-8; %W/m^2*K^4
m = 6*(8.4e-5 + 0.00229112749); %panel (datasheet) + coverglass (density: https://www.tedpella.com/histo_html/coverslip-info.htm)
cap = (310 + 800); %panel (germanium) + coverglass (http://www.markoptics.com/files/Schott%20AF%2032%20eco%20PCP.pdf)
To = 300; %kelvin
period = 5746;
hour = 1:period+1;
xtemp = zeros(1,period+1);
xtemp(1) = To;
% xxtemp = zeros(1,period+1);
% xxtemp(1) = To;
% ytemp = zeros(1,period+1);
% ytemp(1) = To;
% yytemp = zeros(1,period+1);
% yytemp(1) = To;
% xflux = zeros(1,period+1);
% xxflux = zeros(1,period+1);
% yflux = zeros(1,period+1);
% yyflux = zeros(1,period+1);

fluxinx = zeros(1,period+1);
fluxoutx = zeros(1,period+1);

for i = 1:period
    fluxinx(i+1) = qsx(i+1)*alpha*(Asx(i)+fa*At*cos(beta)/pi^2) + qe*epsilon*At/pi;
    fluxoutx(i+1) = -bolt*epsilon*At*xtemp(i)^4;
%     fluxinxx(i+1) = qsxx(i+1)*alpha*(Asxx(i)+fa*At*cos(beta)/pi^2) + qe*epsilon*At/pi;
%     fluxoutxx(i+1) = -bolt*epsilon*At*xxtemp(i)^4;
%     fluxiny(i+1) = qsy(i+1)*alpha*(Asy(i)+fa*At*cos(beta)/pi^2) + qe*epsilon*At/pi;
%     fluxouty(i+1) = -bolt*epsilon*At*ytemp(i)^4;
%     fluxinyy(i+1) = qsyy(i+1)*alpha*(Asyy(i)+fa*At*cos(beta)/pi^2) + qe*epsilon*At/pi;
%     fluxoutyy(i+1) = -bolt*epsilon*At*yytemp(i)^4;
%     
    xflux(i+1) = fluxinx(i+1) + fluxoutx(i+1);
%     xxflux(i+1) = fluxinxx(i+1) + fluxoutxx(i+1);
%     yflux(i+1) = fluxiny(i+1) + fluxouty(i+1);
%     yyflux(i+1) = fluxinyy(i+1) + fluxoutyy(i+1);
%     
    xtemp(i+1) = xflux(i+1) / (m*cap) + xtemp(i);
%     xxtemp(i+1) = xxflux(i+1) / (m*cap) + xxtemp(i);
%     ytemp(i+1) = yflux(i+1) / (m*cap) + ytemp(i);
%     yytemp(i+1) = yyflux(i+1) / (m*cap) + yytemp(i);
end

 
% %K -> C
% for i= 1:period+1
%     xtemp(i) = xtemp(i)-273.15;
% end

days = period+1 / 3600;

plot(hour,xtemp) 
% hold all
% plot(days,xxtemp)
% plot(days,ytemp)
% plot(days,yytemp)
title('Panel temp')
% legend('X Pos', 'X Neg', 'Y Pos', 'Y Neg')
xlabel('Time (s)')
ylabel('Temp (deg K)')

% figure
% plot(hour, fluxin)
% hold all
% plot(hour, fluxout)
% xlabel('Time (s)')
% ylabel('W')