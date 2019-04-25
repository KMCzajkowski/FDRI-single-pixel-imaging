% This is Matlab/Octave function for FS dithering
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
function img=dither(img)
img=img/255;
Nx=size(img,1);
Ny=size(img,2);

img=[zeros(Nx,1) img zeros(Nx,1)];
img=[zeros(1,Ny+2);img;zeros(1,Ny+2)];

for y=2:(Ny+1)
    yimg=img(:,y);
    
    [yimg,qerrv]=octavedith(yimg);
    qerrv=qerrv(2:end-1);
    
    img(:,y)=yimg;
    img(1:Nx,y+1)=img(1:Nx,y+1)+qerrv*3/16;
    img(2:Nx+1,y+1)=img(2:Nx+1,y+1)+qerrv*5/16;
    img(3:Nx+2,y+1)=img(3:Nx+2,y+1)+qerrv*1/16;
end
img=img(2:(Nx+1),2:(Ny+1));
