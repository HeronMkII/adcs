timee = length(xPosFlux) -1;
qsx = xPosFlux; %W/m^2
qsxx = xNegFlux;
qsy = yPosFlux;
qsyy = yNegFlux;
alpha = 0.92;
epsilon = 0.85;
At = solarPanelArea; %m^2
Asx = zeros(1,timee+1); %m^2
Asxx = zeros(1,timee+1);
Asy = zeros(1,timee+1);
Asyy = zeros(1,timee+1);

for i=1:timee+1
    Asx(i) = At*qsx(i)/1395;
    Asxx(i) = At*qsxx(i)/1395;
    Asy(i) = At*qsy(i)/1395;
    Asyy(i) = At*qsyy(i)/1395;

end

fa = 0.227; %????
beta = 20.336*pi/180; %rad
qe = 218; %?????
bolt = 5.67e-8; %W/m^2*K^4
m = 6*(8.4e-5 + 0.00229112749); %panel (datasheet) + coverglass (density: https://www.tedpella.com/histo_html/coverslip-info.htm)
cap = (310 + 800); %panel (germanium) + coverglass (http://www.markoptics.com/files/Schott%20AF%2032%20eco%20PCP.pdf)
To = 300; %kelvin
period = 5746*5;
hour = 1:period+1;
xtemp = zeros(1,period+1);
xtemp(1) = To;
xxtemp = zeros(1,period+1);
xxtemp(1) = To;
ytemp = zeros(1,period+1);
ytemp(1) = To;
yytemp = zeros(1,period+1);
yytemp(1) = To;

tpflux = zeros(1,period+1);

for i = 1:period
    p = [bolt*epsilon*At/(m*cap) 0 0 1 -((qsx(i+1)*alpha*(Asx(i+1) + fa*At*cos(beta)/pi^2) +qe*epsilon*At/pi)/(m*cap) + xtemp(i))];
    xtemp(i+1) = min(roots(p));
    p = [bolt*epsilon*At/(m*cap) 0 0 1 -((qsxx(i+1)*alpha*(Asxx(i+1) + fa*At*cos(beta)/pi^2) +qe*epsilon*At/pi)/(m*cap) + xxtemp(i))];
    xxtemp(i+1) = min(roots(p));
    p = [bolt*epsilon*At/(m*cap) 0 0 1 -((qsy(i+1)*alpha*(Asy(i+1) + fa*At*cos(beta)/pi^2) +qe*epsilon*At/pi)/(m*cap) + ytemp(i))];
    ytemp(i+1) = min(roots(p));
    p = [bolt*epsilon*At/(m*cap) 0 0 1 -((qsyy(i+1)*alpha*(Asyy(i+1) + fa*At*cos(beta)/pi^2) +qe*epsilon*At/pi)/(m*cap) + yytemp(i))];
    yytemp(i+1) = min(roots(p));
    tpflux(i+1) = (xtemp(i+1)-xtemp(i) + xxtemp(i+1)-xxtemp(i) + ytemp(i+1)-ytemp(i) + yytemp(i+1)-yytemp(i)) / (m*cap);
end

days = period.*(1/3600);

powerin = sum(tpflux)

figure
plot(hour,xtemp)
hold on
plot(hour,xxtemp)
plot(hour,ytemp)
plot(hour,yytemp)
title('Panel temp')
xlabel('Time (s)')
ylabel('Temp (deg K)')
legend('X pos','X neg', 'Y pos', 'Y neg')
hold off
figure
plot(hour, tpflux)
title('Heat in')
ylabel('Joules')
xlabel('Time (s)')
