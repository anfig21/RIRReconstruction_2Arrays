function DOA = dirDOA_CS(P,Nnorm,H,uk)
%DOA = dirDOA_CS(Data,Direct,Dict,plotFlag) Applies Compressive
%Sensing to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - P     : frequency response at the specific frequency bins. Nf x M
%       - Nnorm : noise norm estimation. Nf x 1
%       - H     : dictionary of plane waves. M x N x Nf
%       - uk    : unit directional vector. 3 x N
%   Output:
%       - DOA   : DOA estimation via CS. 3 x Nf
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
if nargin < 3, error('dirDOA_CS Error: Not enough input parameters.'), end

%% MAIN CODE
Nf = size(P,1);
N = size(H,2);

DOA = nan(3,Nf);
NoiseMargin = 10;           % dB

c = waitbar(0,'Loading...0\%','Name','dirDOA_CS: CVX across frequencies...');
for ii = 1:Nf
    epsilon = 10^(NoiseMargin/20)*Nnorm(ii);
    Hii = squeeze(H(:,:,ii));
    pii = P(ii,:).';
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
        variable x(N) complex;
        minimize norm(x,1);
        subject to
            norm((Hii*x-pii),2) <= epsilon;
    cvx_end
    
    [~,Idx] = max(abs(x));
    
    DOA(:,ii) = -uk(:,Idx);
    
    waitbar(ii/Nf,c,strcat("Loading... ",string(round(100*ii/Nf,2)),"\,\%"));
end
delete(c)

disp('Direct sound: DOA - CS... OK')

end

