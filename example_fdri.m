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
function example_fdri
if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
end

close all
%----------------------------------
fprintf('Fourier Domain Regularized Inversion (FDRI) example:\n')
N=256;  % all images will be resized to NxN pixels
k=1966; % number of basis functions (the compression ratio is equal to k/N^2)
mi=0.5; % FDRI parameters mi,eps (Eq. 11)
ep=1e-5;
method=0; % =0: FDRI is calculated with pinv() function; =1: with inv(); =2 with svd (see fdri.m)
binarize=true;% use binarized or continuous DCT functions for the measurement matrix
random_basis_selection=false; % random or deterministic selection of DCT bases (See SelectionMatrix.m)
imgdb_path='images/'; % path to an image database with images 1.gif, 2.gif,...
tstimg_path='tst_images/'; % path to test images
tst_images={'lena512.bmp','bird512.jpg','fox512.gif','FUWchart512.jpg'}; % path to the test image

fprintf('Resolution: [%d x %d]\n',N,N);
fprintf('Compression ratio: %.3g%%\n',k/N^2*100);
fprintf('FDRI parameters: mi=%g, eps=%g\n',mi,ep)
txt1={'real-valued','binarized'};
txt2={'deterministically','randomly'};
fprintf('Sampling functions: subset of DCT bases, %s, selected %s\n\n', txt1{1+binarize},txt2{1+random_basis_selection})
dim=[N,N];
fprintf('\nI. PREPARATION STAGE\n');
fprintf('1. Calculate the average DCT spectrum using an image database...\n');
[AvgDCT]=AvgDCTSpectrum(dim,imgdb_path);

fprintf('2. Select %d DCT bases based on the frequency of their occurances in the image databse...\n',k)
[SM]=SelectionMatrix(k,AvgDCT,random_basis_selection);

fprintf('3. Prepare the %dx%d measurement matrix consisting of the selected %s DCT functions...\n',k,prod(dim),txt1{1+binarize})


[M]=MeasurementMatrix(SM,binarize);


fprintf('4. Calculate the generalized inverse matrix with FDRI (takes some time)...\n')
[P]=fdri(M,dim(2),dim(1),mi,ep,method);
M=single(M); % convert the measurement and reconstruction matrices
P=single(P); % to single precision - to save memory and increase reconstruction speed

fprintf('\nII. COMPRESSIVE MEASUREMENTS\n');

for testnr=1:length(tst_images)
    fprintf('1. Prepare scene %s...\n',tst_images{testnr});
    xorig=imresize(double(imread([tstimg_path,tst_images{testnr}])),dim);
    x=single(xorig(:));
    figure(testnr)
    fprintf('2. Take the compressive measurement %d...\n',testnr);
    y=M*x; % Eq. (1)
    fprintf('3. Reconstruct image %d with FDRI...',testnr);
    tic;
    x0=P*y;
    t=toc;
    x0=double(reshape(x0,dim));
    PSNR=psnr(x0,xorig,max(xorig(:)));
    fprintf('(time=%.4gs, psnr=%.4g dB)\n\n',t,PSNR);
    
    subplot(2,3,1)
    imagesc(xorig);colormap(gray);title('Reference image');axis image
    subplot(2,2,2)
    imagesc(double(SM));title('Selected DCT basis')
    subplot(2,2,3)
    plot(y,'.r')
    title(sprintf('Measured data (compr. ratio= %.3g%%)',100*k/N^2))
    subplot(2,2,4)
    x0=double(x0);
    imagesc(x0);colormap(gray);
    title(sprintf('Reconstructed image (FDRI,mi=%.3g,PSNR=%.4g dB)',mi,PSNR));axis image
end
fprintf('\nPlot some examples of sampling functions included in the measurement matrix...\n');
PlotSamplingFunctionExamples(M,dim);
end
