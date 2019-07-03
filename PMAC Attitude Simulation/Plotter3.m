function [] = Plotter3(fignum, t_axis, ary1, ary1_str,unit1, ary2, ary2_str,unit2)
%PLOTTER3 plots a figure with 2 columns for 2 varibles, with 3 rows
%representing each axis
%   Detailed explanation goes here

if isempty(unit1)==0, unit1 = [' (',unit1,')']; end
if isempty(unit2)==0, unit2 = [' (',unit2,')']; end

figure(fignum)
ax = 1;
xlimit = [0, round(t_axis(end),2)];
for indx1 = 1:2:6
    subplot(3,2,indx1)
    plot(t_axis,ary1(ax,:));
    ylabel(sprintf('%s_%d%s',ary1_str,ax,unit1));
    xlabel('Orbit');
    grid on; grid minor;
    ax = ax + 1;
    xlim(xlimit);
end

ax = 1;
for indx2 = 2:2:6
    subplot(3,2,indx2)
    plot(t_axis,ary2(ax,:));
    ylabel(sprintf('%s_%d%s',ary2_str,ax,unit2));
    xlabel('Orbit');
    grid on; grid minor;
    ax = ax + 1;
    xlim(xlimit)
end
end

