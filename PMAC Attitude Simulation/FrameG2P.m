function [C_PG] = FrameG2P(i, o, w)
%FRAMEG2P Frame rotation, Earth inertial frame to Perifocal frame
%   based on orbital parameters (in degrees)
    ci=cosd(i); co=cosd(o); cw=cosd(w);
    si=sind(i); so=sind(o); sw=sind(w);
    
    C_PG = [cw*co-sw*ci*so, cw*so+sw*ci*co, sw*si;
        -sw*co-cw*ci*so, -sw*so+cw*ci*co, cw*si;
        si*so, -si*co, ci];
end

