function [AvgDCT]=AvgDCTSpectrum(dim,img_path,fdtransform)
% Calculate the magnitude of the DCT spectrum averaged over an image database
% img_path is the path to images and images are called 'n.gif' with n=1..49
% Images are first resized to size dim
% fdtransform may be used for other transforms than the DCT. 
% In this case, fdtransform should be the transform to be used
%
% This code is part of the FDRI package
% https://www.igf.fuw.edu.pl/fdri
% Copyright (C) 2018 K. M. Czajkowski, A. Pastuszczak and R. Kotyński
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
if nargin<3
  fdtransform=@(v)dct2(v);
end

try
  AvgDCT=0;
  N=0;
  while(N<=1000)
      x=imresize(double(imread([img_path,num2str(N+1),'.gif'])),dim);
      AvgDCT=AvgDCT+abs(fdtransform(x));
      N=N+1;
  end

  
catch
  if N>0
    fprintf('(Read %d database images 1.gif, 2.gif...)\n',N)
  else
    fprintf('WARNING: cannot read the image database. Please make sure to put a set images (1.gif, 2.gif...) in the %s folder\n',img_path)
    pause (1);
    fprintf('\n... continuing with a simple model for the average image spectrum...\n')
    [x,y]=meshgrid(linspace(1e-3,1,dim(2)),linspace(1e-3,1,dim(1)));
    AvgDCT=1./(x+y);
    AvgDCT=AvgDCT/AvgDCT(1);
  end
end
end
