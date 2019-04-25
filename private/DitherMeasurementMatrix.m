% This is a sample Matlab/Octave program which reconstructs an image measured compressively
% with binarized or continuous DCT functions using the FDRI method.
%
% Krzysztof M. Czajkowski, Anna Pastuszczak, and Rafał Kotyński
% "Real-time single-pixel video imaging with Fourier domain regularization,"
% Optics Express, vol. 26(16), pp. 20009-20022, (2018).
% http://dx.doi.org/10.1364/OE.26.020009
%
%
% This code is part of the FDRI package
% https://www.igf.fuw.edu.pl/fdri
% Copyright (C) 2018/2019 K. M. Czajkowski, A. Pastuszczak and R. Kotyński
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
function M=DitherMeasurementMatrix(M1,dim1,dim2)
if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
    addpath additional    
end
N=size(M1,1);
M= false(N,dim2(1)*dim2(2));
parfor i=1:N
    if ~mod(i,500) || i==1
        disp(['number of binarized function: ',num2str(i),'/',num2str(N)]);
    end
    h=reshape(M1(i,:),dim1);
    if mean(h(:))~=0
        h=h/mean(h(:))*0.5;
        h(h>1)=1;
    end
    test0=myinterp2(255*h,dim2);
    %tic
    h0=dither(test0);
    %toc
    M(i,:)=logical(h0(:));
end
end

function y=myinterp2(x,dim)
dim0=size(x);
x0=linspace(-1,1,dim0(2));
y0=linspace(-1,1,dim0(1));
[X0,Y0]=meshgrid(x0,y0);

x0=linspace(-1,1,dim(2));
y0=linspace(-1,1,dim(1));
[X,Y]=meshgrid(x0,y0);

y=interp2(X0,Y0,x,X,Y);
end