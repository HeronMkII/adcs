function [Across] = skew(A)
%Returns the skew symmetric matrix of a vector in R3
%for use in the multiplication form of the cross product

Across = [0 -A(3) A(2);...
          A(3) 0 -A(1);...
          -A(2) A(1) 0];

end

