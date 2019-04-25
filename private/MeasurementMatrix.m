function [M]=MeasurementMatrix(SM,binarize,invtransform)
% Return the measurement matrix M with rows containing DCT basis selected 
% by the elements of selection matrix SM
% invtransform may be used for other transforms than the DCT. 
% In this case, invtransform should be the inverse transform of the transform used
%
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
dim=size(SM);
m_rows=sum(SM(:));
m_cols=numel(SM);
M=zeros([m_rows,m_cols]);
P=1:m_cols;
P=P(SM);
if nargin<3
  invtransform=@(v)idct2(v);
end
if nargin<2
    binarize=false;
end


if exist ('OCTAVE_VERSION', 'builtin')
    pkg load image
    pkg load signal
    more off
end
% M=dim(1);
% N=dim(2);
% m=1:M;
% n=1:N;
% [n,m]=meshgrid(n,m);
% n=n(:);
% m=m(:);
% put the respective DCT basis in every row of M
  for r=1:m_rows  
    Eye=zeros(dim); 
    Eye(P(r))=1.;
    M(r,:)=reshape(invtransform(Eye),[1,m_cols]);
    %M(r,:)=cos(pi*(2*m+1)*p/(2*M)).*cos(pi*(2*n+1)*q/(2*N));
    if binarize
        M(r,:)=double(M(r,:)>=mean(M(r,:)));
    end
  end
end
