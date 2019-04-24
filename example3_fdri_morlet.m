% This is a sample Matlab/Octave program which reconstructs an image measured compressively
% using the FDRI method. The measurement matrix contains (binarized or real-valued)
% functions obtained as convolutions of Morlet wavelets with white noise
%
% [1]. Krzysztof M. Czajkowski, Anna Pastuszczak and Rafał Kotyński
% "Real-time single-pixel video imaging with Fourier domain regularization,"
% Optics Express, vol. 26(16), pp. 20009-20022, (2018).
% http://dx.doi.org/10.1364/OE.26.020009
% [2]. Krzysztof M. Czajkowski, Anna Pastuszczak and Rafał Kotyński,
% "Single-pixel imaging with Morlet wavelet correlated random patterns,"
% Scientific Reports, vol. 8, art. 466, , (2018).
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function example3_fdri_morlet
if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
end

close all
%----------------------------------
fprintf('Fourier Domain Regularized Inversion (FDRI) example 3:\n')
N=256;  % all images will be resized to NxN pixels
k=fix(N*N*0.03); % number of basis functions (the compression ratio is equal to k/N^2)
mi=0.5; % FDRI parameters mi,eps (Eq. 11)
ep=1e-5;
method=0; % =0: FDRI is calculated with pinv() function; =1: with inv(); =2 with svd (see fdri.m)
tstimg_path='tst_images/'; % path to test images
tst_images={'lena512.bmp','bird512.jpg','fox512.gif','FUWchart512.jpg'}; % path to the test image

binarize=true; % true for the measurement matrix with values {0,1}, false for original WH functions

fprintf('Resolution: [%d x %d]\n',N,N);
fprintf('Compression ratio: %.3g%%\n',k/N^2*100);
fprintf('FDRI parameters: mi=%g, eps=%g\n',mi,ep)
txt1={'real-valued','binarized (0,1)'};
fprintf('Sampling functions: subset of Morlet-wavelets convolved with white noise, %s\n\n', txt1{1+binarize})
dim=[N,N];
fprintf('\nI. PREPARATION STAGE\n');
fprintf('1. Prepare the %dx%d measurement matrix consisting of the selected %s Morlet-wavelets convolved with white noise (takes some time)...\n',k,prod(dim),txt1{1+binarize})

[M,sigma,np,theta]=MorletBasedMeasurementMatrix(k,dim,binarize);

fprintf('2. Calculate the generalized inverse matrix with FDRI (takes some time)...\n')
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
    fprintf('(time=%.3fs, psnr=%.4g dB)\n\n',t,PSNR);
    
    subplot(2,3,1)
    imagesc(xorig);colormap(gray);title('Reference image');axis image
    subplot(2,2,2)
    plot(sigma(2:end)*100,np(2:end),'.k');
    ylabel('Number of peaks n_p')
    xlabel(sprintf('\\sigma (%%) = 100%% x \\sigma (px) x 3 / %d (px)',N))
    title('Selected Morlet wavelet coefficients')
    axis ij tight
    
    subplot(2,2,3)
    plot(y,'.r')
    title(sprintf('Measured data (compr. ratio= %.3g%%)',100*k/N^2))
    subplot(2,2,4)
    x0=double(x0);
    imagesc(x0);colormap(gray);axis image
    title(sprintf('Reconstructed image (FDRI,mi=%5.3g,PSNR=%4.1f dB)',mi,PSNR));axis image
end
fprintf('\n\nPlot some examples of sampling functions included in the measurement matrix...\n');
if exist ('OCTAVE_VERSION', 'builtin')
    PlotSamplingFunctionExamples(M,dim);
else
    p=PlotSamplingFunctionExamples(M,dim,false);
    for i=1:9        
        subplot(3,3,i)
        title(sprintf('Sample (\\theta=%d^{\\circ}, \\sigma=%d%%, n_p=%.1f)',fix(theta(p(i))*180/pi),fix(100*sigma(p(i))),np(p(i)) ));
    end
    
end
end