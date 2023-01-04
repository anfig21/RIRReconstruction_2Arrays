function DOA = dirDOA_RLS(P,H,uk)
%DOA = dirDOA_RLS(H,Dict) Applies Regularised Least Squares solution with
%Tikhonov regularisation and L-Curve method to the Direct Sound Field to
%obtain the DOA of the source.
%   Input:
%       - P     : frequency response at the specific frequency bins. Nf x M
%       - H     : dictionary of plane waves. M x N x Nf
%       - uk    : unit directional vector. 3 x N
%   Output:
%       - DOA   : DOA estimation via RLS. 3 x Nf
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
if nargin < 3, error('dirDOA_RLS Error: Not enough input parameters.'), end

%% MAIN CODE
Nf = size(P,1);

DOA = nan(3,Nf);
for ii = 1:Nf
    [x,~] = reguLeastSquares(squeeze(H(:,:,ii)),P(ii,:).');
    
    [~,Idx] = max(abs(x));
    DOA(:,ii) = -uk(:,Idx);
end

disp('Direct sound: DOA - RLS... OK')

end

