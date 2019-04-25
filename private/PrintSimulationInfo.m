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
function PrintSimulationInfo(k,N,mi,ep,binarize,random_basis_selection)
fprintf('Resolution: [%d x %d]\n',N,N);
fprintf('Compression ratio: %.3g%%\n',k/N^2*100);
fprintf('FDRI parameters: mi=%g, eps=%g\n',mi,ep)
txt1={'real-valued','binarized'};
txt2={'deterministically','randomly'};
fprintf('Sampling functions: subset of DCT bases, %s, selected %s\n\n', txt1{1+binarize},txt2{1+random_basis_selection})
end