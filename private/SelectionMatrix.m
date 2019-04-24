
function [SM]=SelectionMatrix(k,AvgDCT,random_basis_selection)
% Create a (logical) selection matrix SM for the DCT (or other) basis
% with k true elements at the selected DCT basis
% If random_basis_selection=true, selection is random with with the Bernoulli 
% probabilities proportional to AvgDCT 
% If random_basis_selection=false, selection is deterministic and is based on the largest
% average DCT components 
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
  random_basis_selection=true;
end

R=AvgDCT;
if random_basis_selection
  R=rand(size(AvgDCT)).*R;
end
% R=AvgDCT; % instead we could select the largest coefficients deterministically
R(1)=AvgDCT(1); % make sure to include the zeroth frequency
[~,I]=sort(R(:));
SM=zeros(size(AvgDCT),'logical');
SM(I(end-k+1:end))=true;           
end
