function  [yt,yd] = hp(y,lambda)

% Sparse implentation
% Sparsity factor = (5T-6)/(T^2) where T = size(y)
% Picture of the Matrix
%
% HP Detrend time series with Hodrick-Prescott method.
% finds the series yt(:,i) that minimizes:
%
% sum(y(t,i)-yt(t,i))^2 + lambda sum [(yt(t+1,i)-yt(t,i))- ..
% t=1:T                         t=2:T-1   -(yt(t,i)-yt(t-1,i))]^2
%
% for each column i=1,...k in y.
% A larger lambda results in a smoother trend series. 
% For quarterly data Hodrick and Prescott(1980) use phi=1600.
%
% Also returns the series with no trend: yd(:,i)=y(:,i)-yt(:,i)
%                                       i=1,...k
%
% References:
%    
%    Hodrick, Robert J. and Edward C. Prescott,"Post-War
%    U.S. Business Cycles: An Empirical Investigation''
%    Journal of Money, Credit and Banking.
%    Prescott, Edward C.,"Theory Ahead of Business Cycle 
%    Measurement'' QUARTERLY REVIEW, Federal Reserve Bank 
%    of Minneapolis, Fall 1986.
%
% Jesus Fernandez-Villaverde, 4-26-01

    T = length(y);
	K = spdiags([ones(T-2,1) -2*ones(T-2,1) ones(T-2,1)],0:2,T-2,T);
    A = sparse(eye(T)+lambda*(K'*K));
    %spy(A)
 	yt = A\y;
    yd = y-yt;
