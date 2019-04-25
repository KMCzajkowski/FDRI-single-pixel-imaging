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
function ShowReconstruction(testnr,PSNR,t,xorig,SM,x0,y,k,N,mi)
    figure(testnr)
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