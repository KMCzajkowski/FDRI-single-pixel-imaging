% This is a sample Matlab/Octave program which reconstructs an image measured compressively
% using the FDRI method. The measurement matrix contains (binarized or real-valued)
% functions obtained as convolutions of Morlet wavelets with white noise
%
% [1] Krzysztof M. Czajkowski, Anna Pastuszczak and Rafał Kotyński,
% "Single-pixel imaging with Morlet wavelet correlated random patterns,"
% Scientific Reports, vol. 8, art. 466, (2018).
% http://dx.doi.org/10.1038/s41598-017-18968-6
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>

function [M,sigma,np,theta]=MorletBasedMeasurementMatrix(k,dim,binarize)
if nargin<3
    binarize=false;
end

theta= pi.*(rand(k,1)-0.5); % generate random parameters (theta,np,sigma)
sigma=abs(0.33+.165*randn(k,1));% with pdf taken from ref [1]
np=abs(6+4*randn(k,1));

% some alternative parameters:
%npeaks=10;np=rand(k,1)*npeaks;
%sigmin=1/6;sigmax=2.5/6;sigma=sigmin+(sigmax-sigmin)*rand(k,1);% some
%alternative parameters

np(np<.01)=.01;
sigma(sigma<0.01)=0.01;

M=zeros(k,prod(dim));
if binarize
    M(1,:)=1;
else
    M(1,:)=1/sqrt(prod(dim)); % make sure to include a constant sampling function
end

x = linspace(- 1.5, 1.5, dim(2));
y = linspace(- 1.5, 1.5, dim(1));
[X, Y] = meshgrid(x, y);

for i=2:k
    % Find a centered Morlet h wavelet with parameters (np(i),sigma(i),theta(i))
    h1 = exp((0.5i*pi*np(i)/sigma(i)) * (X*cos(theta(i))+Y*sin(theta(i))));
    h2 = exp((-X.^2-Y.^2)/(2*sigma(i)^2));
    h=h1.*(h2-sum(h1(:).*h2(:))/sum(h1(:)));
    % Convolve h with white noise
    h=real(ifft2(fft2(h).*exp(2i*pi*rand(size(h)))));
    if binarize
        M(i,:)=reshape(real(real(h)>mean(real(h(:)))),[1,prod(dim)]);
    else
        M(i,:)=reshape(real(h(:)),[1,prod(dim)]);
        M(i,:)=M(i,:)-mean(M(i,:));
        M(i,:)=M(i,:)/norm(M(i,:));
    end
end
