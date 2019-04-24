function [P]=fdri(M,Nx,Ny,mi,ep,method)
% Calculates the generalised inverse of the measurement matrix using the FDRI method
% Czajkowski et al, Opt. Express 26, 20009, 2018, http://dx.doi.org/10.1364/OE.26.020009
%
% Input parameters:
% M - measurement matrix. Every row includes one Nx x Ny sampling function
% mi - FDRI parameter (defaults to mi=0.5)
% ep - FDRI parameter (defaults to ep=1e-5)
% method =0 (default) calculate P with Eq. (7)
% method =1 calculate P with Eq. (8) using the pinv function
% method =2 calculate P with Eq. (8) using svd
%
% Output parameters:
% P - the generalized inverse matrix (calculated with Eq. (7) or (8))
%
% 
% Krzysztof M. Czajkowski K.M., Anna Pastuszczak A., and Rafał Kotyński
% "Real-time single-pixel video imaging with Fourier domain regularization,"
% Optics Express, vol. 26(16), pp. 20009-20022, (2018).
% http://dx.doi.org/10.1364/OE.26.020009
%
%
% This code is part of the FDRI/SIMPLEX package
% https://www.igf.fuw.edu.pl/fdri
% Copyright (C) 2018 K. M. Czajkowski, A. Pastuszczak and R. Kotyński
%
%
%  GPL LICENSE INFORMATION
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


if nargin<6,method=0;end;
if nargin<5,ep=1e-5;end;
if nargin<4,mi=0.5;end;

k=numel(M)/Nx/Ny;
M=reshape(M,[k,Nx*Ny]);

if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
end

% calculate the diagonal elements of hat(Gamma) according to Eq. (11)
[wx,wy]=meshgrid(2*pi/Nx*[0:Nx/2-1,-Nx/2:-1],2*pi/Ny*[0:Ny/2-1,-Ny/2:-1]);
D=1./sqrt( (1-mi)^2*(sin(wx).^2+sin(wy).^2) +ep +mi^2*(wx.^2+wy.^2)/(2*pi^2));

% helper functions - apply 2D DFT to images stored in rows or columns of a matrix X
row_fft2=@(X)reshape(fft2(reshape(X.',[Ny,Nx,size(X,1)])),[Nx*Ny,size(X,1)]).';% size of X is [k, n*n], fft2 is applied to rows
row_ifft2=@(X)reshape(ifft2(reshape(X.',[Ny,Nx,size(X,1)])),[Nx*Ny,size(X,1)]).';% size of X is [k, n*n], ifft2 is applied to rows
col_fft2=@(X)reshape(fft2(reshape(X,[Ny,Nx,size(X,2)])),[Nx*Ny,size(X,2)]);% size of X is [n*n,k], fft2 is applied to columns
col_ifft2=@(X)reshape(ifft2(reshape(X,[Ny,Nx,size(X,2)])),[Nx*Ny,size(X,2)]);% size of X is [n*n,k], ifft2 is applied to columns

% helper functions - apply 2D linear filtering to a matrix X
FILT_R=@(X)row_fft2(row_ifft2(X)*spdiags(D(:),0,Nx*Ny,Nx*Ny)); % X*F'*D*F
FILT_L=@(X)col_fft2(spdiags(D(:),0,Nx*Ny,Nx*Ny)*col_ifft2(X)); % F*D*F'*X

% Now calculate the generalized inverse matrix P
a=real(FILT_R(M));
switch method
    case 1 % calculate the inversion matrix with Eq. (8)
        P=real(FILT_L(pinv(a)));
   
    case 2 %use svd to calculate the pseudoinverse, and then use Eq. (8):
        [U, S, V] = svd (a,'econ'); %a = U*S*V'
        tol=1e-5;
        S=diag(S);
        S(abs(S)<max(size(S)) * max(abs(S(:)))*eps *tol)=0;
        INV_S=1./S;
        INV_S(isnan(INV_S))=0;
        P=real(FILT_L(V*diag(INV_S)*U'));
        
    otherwise % calculate the inversion matrix with Eq. (7) - default
        P=real(FILT_L(a'*(inv(a*a'))));
end


  
