function inv_A = inv_blk_diag(A)
%INV_BLK_DIAG   Computes the inverse of a block diagonal matrix.
%
%   INV_A = INV_BLK_DIAG(A) returns the inverse of a block diagonal matrix.
%   A is an N-by-1 cell array of diagonal block matrices.
%
% Author: Minjie Fan, 2015

inv_A = cellfun(@inv, A, 'UniformOutput', false);

end