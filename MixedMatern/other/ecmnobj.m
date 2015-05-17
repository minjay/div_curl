function Objective = ecmnobj(Data, Mean, Covar, CholCovar)
%ECMNOBJ Evaluate negative log-likelihood function for ECMNMLE.
%	Observed negative log-likelihood function based on current maximum
%	likelihood parameter estimates.
%
%   Objective = ecmnobj(Data, Mean, Covar);
%   Objective = ecmnobj(Data, Mean, Covar, CholCovar);
%
% Inputs:
%   Data - NUMSAMPLES x NUMSERIES matrix of observed multivariate normal data
%		with NaNs to represent missing values.
%   Mean - NUMSERIES x 1 column vector with current mean estimate of Data.
%   Covar - NUMSERIES x NUMSERIES matrix with current covariance estimate of
%		Data.
%
% Optional Inputs:
%   CholCovar - Cholesky decomposition of covariance matrix, i.e., chol(Covar).
%
% Outputs:
%   Objective - Value of the observed negative log-likelihood function over
%		the Data given current estimates for the Mean and Covar of the data.
%
% See also ECMNMLE, ECMNSTD, ECMNFISH, ECMNHESS.

%	Author(s): R.Taylor, 4-21-2005
%	Copyright 2005 The MathWorks, Inc.
%	$Revision: 1.1.6.5 $   $Date: 2012/08/21 00:08:26 $

% Step 1 - check arguments

if nargin < 3
	error(message('finance:ecmnobj:MissingInputArg'));
else
	if isempty(Data)
		error(message('finance:ecmnobj:EmptyInputData'));
	end
	if isempty(Mean)
		error(message('finance:ecmnobj:EmptyInputMean'));
	end
	if isempty(Covar)
		error(message('finance:ecmnobj:EmptyInputCovar'));
	end
	
	[NumSamples, NumSeries] = size(Data);

	Mean = Mean(:);
	if ~all(size(Mean) == [NumSeries, 1])
		error(message('finance:ecmnobj:IncompatibleMean'));
	end
	if ~all(size(Covar) == [NumSeries, NumSeries]) 
		error(message('finance:ecmnobj:IncompatibleCovar'));
	end
end
if nargin < 4 || isempty(CholCovar) || ~all(size(CholCovar) == size(Covar))
	[CholCovar, CholState] = chol(Covar);
	if CholState > 0
		error(message('finance:ecmnobj:NonPosDefCovar'));
	end
end

LogDetCovar = 2.0 * sum(log(diag(CholCovar)));

% Step 2 - initialization

LogTwoPi = log(2.0 * pi);

% Step 3 - main loop over data records

Objective = 0.0;
Count = 0;

for k = 1:NumSamples

	P = ~isnan(Data(k,:));
	Available = sum(P);

	if Available > 0				% skip over empty records
        Count = Count + 1;

		Objective = Objective + 0.5 * Available * LogTwoPi;
	
        if Available < NumSeries
			[SubChol, CholState] = chol(Covar(P,P));
			
			if CholState > 0
				error(message('finance:ecmnobj:NonPosDefSubCovar'));
			end

			SubResid = SubChol' \ (Data(k,P)' - Mean(P));
			
			Objective = Objective + 0.5 * SubResid' * SubResid;
			Objective = Objective + sum(log(diag(SubChol)));
		else
			Resid = CholCovar' \ (Data(k,:)' - Mean);

			Objective = Objective + 0.5 * Resid' * Resid;
			Objective = Objective + 0.5 * LogDetCovar;
		end
    end
end

if Count < 1
    error(message('finance:ecmnobj:AllNaNData'));
end
