function H_A = H_blk_diag(A)
%H_BLK_DIAG   Computes the conjugate transpose of a block diagonal matrix.
%
%   H_A = H_BLK_DIAG(A) returns the conjugate transpose of a block diagonal
%   matrix. A is an N-by-1 cell array of diagonal block matrices.
%
% Author: Minjie Fan, 2015

H_A = cellfun(@ctranspose, A, 'UniformOutput', false);

end