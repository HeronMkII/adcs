function e = eclipsef(t,ec,time,orbits,edur,period)
e=false;
found=false;
% c=1;
% while (c<time && found==false)
%     if ec(c) == 1
%         if t==c
%             found=true;
%             e=true;
%         end
%     end
%     c=c+1;
% end

for k=1:orbits
%     for each orbit
    a = k*period-edur; %beginning of eclipse
    b = k*period; %end of eclipse
    if (a<t)
        if (t<b)
            e=true;
        end
    end
    if (e)
        break
    end
end