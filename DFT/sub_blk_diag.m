function C = sub_blk_diag(A, B)
%SUB_BLK_DIAG   Subtraction of two block diagonal matrices.
%
%   C = SUB_BLK_DIAG(A, B) returns the subtraction of two block diagonal
%   matrices A and B. They are N-by-1 cell arrays of diagonal block
%   matrices.
%
% Author: Minjie Fan, 2015

C = cellfun(@(X, Y) X-Y, A, B, 'UniformOutput', false);

end
