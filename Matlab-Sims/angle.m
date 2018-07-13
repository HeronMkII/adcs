function [theta] = angle(A,B)
%Given 2 column vectors in R3, finds the angle between them in degrees

theta = acosd(A'*B/(norm(A)*norm(B)));

end

