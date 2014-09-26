% bernoull() - Computes Bernoulli distribution of x for 
% 		"natural parameter" eta. The mean m of a 
%		Bernoulli distributions relates to eta as,
% 		m = exp(eta)/(1+exp(eta));
%
% Usage:
%   >> [p]=bernoull(x,eta);
%
% Inputs:
%   x		- data
%   eta		- distribution parameter
%
% Outputs:
%   p		- probability
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)
%          and Paul Sajda (ps629@columbia,edu 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson, Lucas Parra and Paul Sajda
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [p]=bernoull(x,eta);

e = exp(eta);

p = ones(size(e));  

indx = find(~isinf(e));

p(indx) = exp(eta(indx).*x - log(1+e(indx)));
