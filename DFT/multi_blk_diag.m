function C = multi_blk_diag(A, B)
%MULTI_BLK_DIAG   Multiplication of two block diagonal matrices.
%
%   C = MULTI_BLK_DIAG(A, B) returns the multiplication of two block
%   diagonal matrices A and B. They are N-by-1 cell arrays of diagonal
%   block matrices.
%
% Author: Minjie Fan, 2015

C = cellfun(@(X, Y) X*Y, A, B, 'UniformOutput', false);

end