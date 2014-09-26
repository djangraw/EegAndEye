% rocarea() - computes the area under the ROC curve
% 		If no output arguments are specified
%		it will display an ROC curve with the
%		Az and approximate fraction correct.
%
% Usage:
%  >> [Az,tp,fp,fc]=rocarea(p,label);
%
% Inputs:
%   p		- classification output
%   label	- truth labels {0,1}
%
% Outputs:
%   Az		- Area under ROC curve
%   tp		- true positive rate
%   fp		- false positive rate
%   fc		- fraction correct
%
% Authors: Lucas Parra (parra@ccny.cuny.edu, 2004)
%	   with Adam Gerson (reformatted for EEGLAB)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Lucas Parra
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

function [Az,tp,fp,fc]=rocarea(p,label);

[tmp,indx]=sort(-p);

label = label>0;

Np=sum(label==1);
Nn=sum(label==0);

tp=0; pinc=1/Np;
fp=0; finc=1/Nn;
Az=0;

N=Np+Nn;

tp=zeros(N+1,1);
fp=zeros(N+1,1);

for i=1:N
  
  tp(i+1)=tp(i)+label(indx(i))/Np;
  fp(i+1)=fp(i)+(~label(indx(i)))/Nn;
  Az = Az + (~label(indx(i)))*tp(i+1)/Nn;

end;

[m,i]=min(fp-tp);
fc = 1-mean([fp(i), 1-tp(i)]);

if nargout==0
  plot(fp,tp); axis([0 1 0 1]); hold on
  plot([0 1],[1 0],':'); hold off
  xlabel('false positive rate') 
  ylabel('true positive rate') 
  title('ROC Curve'); axis([0 1 0 1]); 
  text(0.4,0.2,sprintf('Az = %.2f',Az))
  text(0.4,0.1,sprintf('fc = %.2f',fc))
  axis square
end
  












