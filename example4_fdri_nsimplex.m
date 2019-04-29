% This is a sample Matlab/Octave program which reconstructs an image measured compressively
% with binarized or continuous DCT functions. The sampling functions are
% distributed over simplex vertices. This enables to reduce impact of
% experimental noise on the reconstruction and remove constant bias from
% the measurements.
%
% Krzysztof M. Czajkowski, Anna Pastuszczak, and Rafał Kotyński
% "Single-pixel imaging with sampling distributed over simplex vertices ,"
% 
%  Optics Letters Vol. 44, Issue 5, pp. 1241-1244 (2019)
%  https://doi.org/10.1364/OL.44.001241 
%
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
function example4_fdri_nsimplex
if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
    addpath additional
end
close all
%----------------------------------
fprintf('Fourier Domain Regularized Inversion (FDRI) example:\n')
N=256;  % all images will be resized to NxN pixels
k=1966; % number of basis functions (the compression ratio is equal to k/N^2)
p=10; % word size (p-before coding, p+1 - after coding)
dim2=[1024 768]; % DMD dimensions
mi=0.5; % FDRI parameters mi,eps (Eq. 11)
ep=1e-5;
method=0; % =0: FDRI is calculated with pinv() function; =1: with inv(); =2 with svd (see fdri.m)
binarize=true;% use binarized or continuous simplex coded functions
savematrix.flag=true;
savematrix.filename='test.bin';
random_basis_selection=false; % random or deterministic selection of DCT bases (See SelectionMatrix.m)
imgdb_path='images/'; % path to an image database with images 1.gif, 2.gif,...
tstimg_path='tst_images/'; % path to test images
tst_images={'lena512.bmp','bird512.jpg','fox512.gif','FUWchart512.jpg'}; % test image names

PrintSimulationInfo(k,N,mi,ep,false,random_basis_selection)

dim1=[N,N];
fprintf('\nI. PREPARATION STAGE\n');
fprintf('1. Calculate the average DCT spectrum using an image database...\n');
[AvgDCT]=AvgDCTSpectrum(dim1,imgdb_path);

fprintf('2. Select %d DCT bases based on the frequency of their occurances in the image databse...\n',k)
[SM]=SelectionMatrix(k,AvgDCT,random_basis_selection);

fprintf('3. Prepare the %dx%d measurement matrix consisting of the selected DCT functions...\n',k,prod(dim1))
[M]=MeasurementMatrix(SM);

fprintf('4. Coding of sampling functions on p-simplex with p=%d \n',p)
[M_enc,Q]=SimplexCoding(M,p);
%M_enc=single(M_enc); % convert the measurement and reconstruction matrices
clear M % remove original sampling matrix to save memory

if binarize
    % The functions are scaled up to DMD size and dithered with error
    % diffusion algorithm. Binarization uses parallel for loops if parallel
    % computing toolbox is available.
    fprintf('5. Dithering of sampling functions and preparing a resized version of the dithered matrix to DMD size %dx%d\n',dim2(1),dim2(2))
    M_enc=DitherMeasurementMatrix(M_enc,dim1,dim2);
    M_enc_small=Downsize(M_enc);
    % Then, FDRI is calculated on a downscaled version of the sampling
    % matrix
    fprintf('6. Calculate the generalized inverse matrix with FDRI (takes some time)...\n')
    [P]=fdri(Q*M_enc_small,dim1(2),dim1(1),mi,ep,method);
    clear M_enc_small
else
    fprintf('5. Calculate the generalized inverse matrix with FDRI (takes some time)...\n')
    [P]=fdri(Q*M_enc,dim1(2),dim1(1),mi,ep,method);
end
P=single(P); % convert to single precision to save memory and increase reconstruction speed
Q=single(Q);
P=P*Q;
clear Q
fprintf('\nII. COMPRESSIVE MEASUREMENTS\n');
%%
for testnr=1:length(tst_images)
    fprintf('1. Prepare scene %s...\n',tst_images{testnr});
    xorig=LoadImageAndResize([tstimg_path,tst_images{testnr}],dim1);
    xorigb=LoadImageAndResize([tstimg_path,tst_images{testnr}],dim2);
    if binarize
        x=single(xorigb(:)); % if binary sampling is used the dimension of
        % sampling function is equal to that of DMD (dim2)
    else
        x=single(xorig(:));  % else it is the same as the original sampling function size
    end
    fprintf('2. Take the compressive measurement %d...\n',testnr);
    y=M_enc*x; % Eq. (1)
    fprintf('3. Reconstruct image %d with FDRI...',testnr);
    tic;
    x0=P*y;
    t=toc;
    x0=double(reshape(x0,dim1));
    PSNR=psnr(x0,xorig,max(xorig(:)));
    ShowReconstruction(testnr,PSNR,t,xorig,SM,x0,y,k,N,mi)
end
fprintf('\nPlot some examples of sampling functions included in the measurement matrix...\n');
if binarize
    PlotSamplingFunctionExamples(M_enc,dim2);
else
    PlotSamplingFunctionExamples(M_enc,dim);
end
if savematrix.flag
    fprintf('\nSaving the measurement matrix...\n');
    SaveSamplingMatrix(M_enc,savematrix.filename)
end
end
