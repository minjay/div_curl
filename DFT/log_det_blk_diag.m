function log_det_A = log_det_blk_diag(A)
%LOG_DET_BLK_DIAG   Computes the log determinant of a block diagonal matrix.
%
%   LOG_DET_A = LOG_DET_BLK_DIAG(A) returns the log determinant of a block
%   diagonal matrix. A is an N-by-1 cell array of diagonal block matrices.
%
% Author: Minjie Fan, 2015

log_det_A = sum(log(cellfun(@det, A)));

end
