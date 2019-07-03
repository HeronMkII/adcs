function [w_cross] = cross_mat(w)
%cross_mat.mat computes the matrix for cross product operations
%   (only for 3-row column vector (3D vectors))
    
    w1 = w(1); w2 = w(2); w3 = w(3);
    
    w_cross = [0 -w3 w2;
                w3 0 -w1;
                -w2 w1 0];

end