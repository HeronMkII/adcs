function [] = PlotterStruct(fignum, var,t_axis,unit)
%PLOTTERSTRUCT plots all the field array inside a struct
%   Detailed explanation goes here
    
    var_ary = fieldnames(var);
    %var_ary = [var_ary(1);var_ary(3:end)]; %removes impulse plot
    numfield = length(var_ary);
    xlimit = [0, round(t_axis(end),2)];

    figure(fignum)
    
    for indx2 = 1:3
        subplot(3,1,indx2)
        hold on
        xlabel('Orbit')
            ylabel(sprintf('torq_%d(%s)',indx2,unit))
            xlim(xlimit)
            grid on; grid minor;
        for indx1 = 1:numfield-1
            ary = getfield(var,string(var_ary(indx1)));
            plot(t_axis,ary(indx2,:)) 
        end
        hold off
        legend(var_ary(1:numfield-1))
    end
    
    
end

