function [fw] = S(w)
%Given the angular velocity of a body, returns the matrix S(w)
%required in the equation dq/dt = 0.5*S(w)*q

fw = [0 -w(1) -w(2) -w(3);...
      w(1) 0 w(3) -w(2);...
      w(2) -w(3) 0 w(1);...
      w(3) w(2) -w(1) 0];


end

