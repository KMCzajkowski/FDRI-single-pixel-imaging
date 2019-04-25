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
function [Msmall]=Downsize(M,Nx,Ny,resix,resiy)
if nargin<2
  resiy=1024;
  resix=768;
  Nx=256;
  Ny=256;
end
if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
end
Msmall=zeros(size(M,1),size(M,2)/(resix/Nx*resiy/Ny));

downsizerow=@(M)myinterp2(conv2(M,ones(resiy/Ny,resix/Nx)/(resix/Nx*resiy/Ny),'same'),[Nx Ny]);
parfor r=1:size(M,1)
    h=downsizerow( reshape(M(r,:),resiy,resix));
    Msmall(r,:)= reshape(h,1,[]);
end
fac=(resix*resiy)/(Nx*Ny);
Msmall=Msmall*fac;
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
