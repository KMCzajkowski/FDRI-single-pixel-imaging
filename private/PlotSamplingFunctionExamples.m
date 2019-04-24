% Plot 9 randomly selected sampling functions from the measurement matrix M
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

function [P]=PlotSamplingFunctionExamples(M,dim,includetitle)
if nargin<3
    includetitle=true;
end

figure();
P=zeros(1,9);
for nr=1:9
    subplot(3,3,nr)
    p=randi(size(M,1));
    imagesc(reshape(real(double((M(p,:)))),dim));axis image
    colorbar();
    if includetitle
        title(sprintf('sampling function %d of %d',p,size(M,1)));
    end
    P(nr)=p;
end
end