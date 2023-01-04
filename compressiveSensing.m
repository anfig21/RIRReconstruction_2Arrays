function x = compressiveSensing(H,p,epsilon)
%x = compressiveSensing(H,N,p,epsilon) Solves the optimisation problem for
%the coefficients x in
%                   min ||p - Hx||  s.t. ||x||_1
%assuming a sparse structure using the cvx toolbox for a given frequency.
%   Input:
%       - H         : sensing matrix. M x N
%       - p         : measured data. M x 1
%       - f         : frequency span. 1 x Nf
%       - epsilon   : estimated variance of the noise. Scalar
%   Output:
%       - x         : estimated coefficients. N x 1
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
if nargin < 3, error('compressiveSensing Error: Not enough input arguments.'), end

[M,N] = size(H);
if M ~= numel(p), error('compressiveSensing Error: Data and sensing matrix dimensions do not agree.'), end

%% MAIN CODE
p = p(:);

cvx_begin quiet
    cvx_precision high
    variable x(N) complex;
    minimize norm(x,1);
    subject to
        norm((H*x-p),2) <= epsilon;
cvx_end

end

