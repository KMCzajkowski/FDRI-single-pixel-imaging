% This is a sample Matlab/Octave program which performs encoding of
% sampling functions on N-simplex
%
% Krzysztof M. Czajkowski, Anna Pastuszczak, and Rafał Kotyński
% "Single-pixel imaging with sampling distributed over simplex vertices ,"
% 
%  Optics Letters Vol. 44, Issue 5, pp. 1241-1244 (2019)
%  https://doi.org/10.1364/OL.44.001241 
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
function [M_enc,D]=SimplexCoding(M,p)
[V,invV]=nt_vortices(p); % find vortives of the n-tetrahedron
M_enc=nt_encoder(V,invV,M); % calculate the self-balanced non-negative measurement matrix
D=nt_decoder(V,size(M,1)); % calculate the decoding matrix
end

function [D]=nt_decoder(V,n)
% V - n-tetrahedron vertices points (in columns)
% n - number of the sampling functions (before coding)
% D - decoding matrix D=kron(I,V)
%
% The modified reconstruction matrix is P1=P*D, 
% where P is the generalized inverse of the measurement matrix
% Parameter n is equal to size(P,2)

p=size(V,1); % word size (p-before coding, p+1 - after coding)
k=fix(n/p)+(mod(n,p)>0); % number of words used in coding 
D=kron(eye(k),V); % block matrix with blocks V
D=D(1:n,:); % truncate the decoding matrix
end

function [M_nt,M]=nt_encoder(V,invV,M)
% M - the measurement matrix_type
% [V,invV]=nt_vortices(p)
% M_nt - enlarged measurement matrix with nonnegative values
% (rows of M contain the sampling functions)

p=size(V,1); % word size (p-before coding, p+1 - after coding)
n=size(M,1); %number of the sampling functions (before coding)
k=fix(n/p)+(mod(n,p)>0); % number of words used in coding (for one column of M)
N=size(M,2); % resolution (rows x columns)
M_nt=zeros(k*(p+1),N);  
if n<k*p
  M=[M;zeros(k*p-n,N)];
end
 
M=reshape(M,p,k*N) ;

M_map=zeros(1,k*N);
dist_map=zeros(1,k*N);

% Create a map of M with the words containing sampling points that should 
% be coded with particular bases
% This is equivalent to finding the most distant point from V
for ii=1:p+1 
  d= sum( abs(V(:,ii) - M).^2 ,1);
  I=(d>dist_map);
  M_map(I)=ii;
  dist_map(I)=d(I);   
end


for ii=1:p+1
  I=(M_map==ii);
  if sum(I(:))
    U=invV(:,:,ii)* reshape(M(logical(kron(I,ones(p,1)))),p,[]);%!!!!kron
    M_nt(logical(kron(I,ones(p+1,1))))=U(:);
  end
end
M_nt=reshape(M_nt,(p+1)*k,N);
M=reshape(M,p*k,N);
end

function [V,invV]=nt_vortices(p)
% Generate a [p x (p+1)] matrix V with columns pointing to the vertices of 
% a p-simplex (p=2 for a triangle, p=3 for a tetrahedron etc) 
% in a p-dimensional space 
% if invV is requested, the set of (p+1) inverse matrices is calculated. 
% The size of invV is [p+1,p,p+1] and invV(:,:,c) is calculated as the inverse
% of V with column c removed and with a row of zeros inserted to invV(:,:,c) 
% at row c

V=zeros(p,p+1);
for ii=1:p
  c=-1/ii;
  s=sqrt(1-c^2);
  V(1:ii-1,1:ii)=V(1:ii-1,1:ii)*s;
  V(ii,1:ii)=c;
  V(ii,ii+1)=1;
end
if nargout>1
  invV=zeros([p+1,p,p+1]);
  insert_empty_row=@(A,r)[A(1:r-1,:);zeros(1,size(A,2));A(r:size(A,1),:)]; %insert an empty row into a matrix
  omit_column=@(A,c)A(:,[1:c-1,c+1:size(A,2)] );
  for ii=1:p+1
    invV(:,:,ii)=insert_empty_row(inv(omit_column(V,ii)),ii);
  end
end
end
