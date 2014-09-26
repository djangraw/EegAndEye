% pop_iirfilt() - interactively filter EEG dataset data using iirfilt()
%
% Usage:
%   >> EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw);
%
% Graphical interface:
%   "Lower edge ..." - [edit box] Lower edge of the frequency pass band (Hz) 
%                 Same as the 'locutoff' command line input.
%   "Higher edge ..." - [edit box] Higher edge of the frequency pass band (Hz) 
%                 Same as the 'hicutoff' command line input.
%   "Notch filter" - [edit box] provide the notch range, i.e. [45 55] 
%                 for 50 Hz). This option overwrites the low and high edge limits
%                 given above. Set the 'locutoff' and 'hicutoff' values to the
%                 values entered as parameters, and set 'revfilt to 1, to swap
%                 from bandpass to notch filtering.
%   "Filter length" - [edit box] Filter lenghth in point (default: see 
%                 >> help pop_iirfilt). Same as 'trans_bw' optional input.
%
% Inputs:
%   EEG       - input dataset
%   locutoff  - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz) {0 -> highpass}
%   trans_bw  - length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt   - [0|1] Reverse filter polarity (from bandpass to notch filter).
%                     Default is 0 (bandpass).
%
% Outputs:
%   EEGOUT   - output dataset
%
% Authors: Maksym Pozdin (mpozdin.ece04@gtalumni.org, IOL/ONRC,2004), 
%          with Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD, La Jolla CA)
%
% See also: iirfilt(), pop_eegfilt(), eegfilt(), eegfiltfft(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: pop_iirfilt.m,v $
% Revision 1.6  2005/10/27 17:09:52  arno
% nothing
%
% Revision 1.5  2005/09/30 17:06:06  arno
% added revfilt etc...
%
% Revision 1.4  2005/09/30 17:04:48  arno
% revision 1.1
%
% Revision 1.1  2005/09/30 16:50:28  arno
% Initial revision
%
% Revision 1.23  2003/12/03 18:31:32  arno
% implementing eegfiltfft
%
% Revision 1.22  2003/09/01 18:14:17  arno
% fixing nargin problem -thanks Petr Janata
%
% Revision 1.21  2003/08/06 00:22:05  arno
% removing debug message
%
% Revision 1.20  2003/08/02 21:29:50  arno
% text
%
% Revision 1.19  2003/07/28 17:39:07  arno
% nothing
%
% Revision 1.18  2003/07/22 17:36:44  arno
% subtract DC for small portion of data
%
% Revision 1.17  2003/07/20 19:19:10  scott
% clarify processing message
%
% Revision 1.16  2003/04/24 22:13:51  arno
% typo
%
% Revision 1.15  2003/04/24 22:13:17  arno
% removing 0 warning
%
% Revision 1.14  2003/04/24 22:08:55  arno
% updating error message
%
% Revision 1.13  2003/04/24 22:06:31  arno
% show message when processing data portion
%
% Revision 1.12  2003/04/24 22:00:16  arno
% debuging boundaries
%
% Revision 1.11  2003/02/17 02:43:43  arno
% reformating text for new functionality in help2html
%
% Revision 1.10  2003/02/16 23:10:43  arno
% adding GUI info
%
% Revision 1.9  2003/01/24 04:03:37  scott
% edits msgs -sm
%
% Revision 1.8  2003/01/24 00:23:35  arno
% debugged revfilt parameter
%
% Revision 1.7  2002/11/15 01:45:53  scott
% can not -> cannot
%
% Revision 1.6  2002/10/16 21:42:25  arno
% default for highcuroff
%
% Revision 1.5  2002/08/12 16:25:02  arno
% inputdlg2
%
% Revision 1.4  2002/08/09 01:43:18  arno
% [6~[6~same
%
% Revision 1.3  2002/08/09 01:42:40  arno
% debugging filter over smal time period
%
% Revision 1.2  2002/08/09 00:41:22  arno
% updating for boundaries
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% Updated 2/26/14 by DJ - added filter outputs b,a

function [EEG, com, b, a] = pop_iirfilt_returnfilter( EEG, locutoff, hicutoff, trans_bw, revfilt);

com = '';
if nargin < 1
	help pop_iirfilt;
	return;
end;	
if isempty(EEG.data)
    disp('pop_iirfilt() error: cannot filter an empty dataset'); return;
end;    
if nargin < 2
	% which set to save
	% -----------------
   	promptstr = { 'Highpass: low edge of the frequency pass band (Hz) (0 -> lowpass)', ...
   				  'Lowpass: high edge of the frequency pass band (Hz) (0 -> highpass)', ...
   				  strvcat('Notch filter the data. Give the notch range, i.e. [45 55] for 50 Hz)', ...
                  '(this option overwrites the low and high edge limits given above)'), ...
                  'Transition BW filter (default: see >> help pop_iirfilt)' };
	inistr       = { '0', '0', '', '' };
   	result       = inputdlg2( promptstr, 'Filter the data -- pop_iirfilt()', 1,  inistr, 'pop_iirfilt');
	if size(result, 1) == 0 return; end;
	locutoff   	 = eval( result{1} );
	hicutoff 	 = eval( result{2} );
	if isempty( result{3} )
		 revfilt = 0;
	else 
        revfilt    = eval( [ '[' result{3} ']' ] );
        locutoff = revfilt(1);
        hicutoff = revfilt(2);
        revfilt = 1;
	end;
	if locutoff == 0 & hicutoff == 0 return; end;
	if isempty( result{4} )
		 trans_bw = [];
	else trans_bw    = eval( result{4} );
	end;
else
    if nargin < 3
        hicutoff = 0;
    end;
    if nargin < 4
        trans_bw = [];
    end;
    if nargin < 5
        revfilt = 0;
    end;
end;
 
options = { EEG.srate, locutoff, hicutoff, EEG.pnts };
if ~isempty( trans_bw )
	options = { options{:} trans_bw };
else 
	options = { options{:} 0 };
end;
if revfilt ~= 0
	options = { options{:} revfilt };
end;

% warning
% -------
if exist('filtfilt') ~= 2
    disp('Warning: cannot find the signal processing toolbox');
    disp('         a simple fft/inverse fft filter will be used');
end;

if EEG.trials == 1 
	if ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
		boundaries = strmatch('boundary', {EEG.event.type});
		if isempty(boundaries)
            [EEG.data,b,a] = iirfilt( EEG.data, options{:}); 
		else
			options{4} = 0;
			disp('pop_iirfilt:finding continuous data boundaries');
			tmplat = cell2mat({EEG.event.latency});
            boundaries = tmplat(boundaries);
            boundaries = [0 round(boundaries-0.5) EEG.pnts];
            try, warning off MATLAB:divideByZero
            catch, end;
			for n=1:length(boundaries)-1
				try
                    fprintf('Processing continuous data (%d:%d)\n',boundaries(n),boundaries(n+1)); 
                    [EEG.data(:,boundaries(n)+1:boundaries(n+1)),b,a] = ...
                        iirfilt(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{:});
				catch
					fprintf('\nFilter error: continuous data portion too narrow (DC removed if highpass only)\n');
                    if locutoff ~= 0 & hicutoff == 0
                        tmprange = [boundaries(n)+1:boundaries(n+1)];
                        EEG.data(:,tmprange) = ...
                            EEG.data(:,tmprange) - repmat(mean(EEG.data(:,tmprange),2), [1 length(tmprange)]);
                    end;
				end;
			end
            try, warning on MATLAB:divideByZero
            catch, end;
		end
	else
        [EEG.data,b,a] = iirfilt( EEG.data, options{:});
	end;
else
    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
    [EEG.data,b,a] = iirfilt( EEG.data, options{:});
	% Note: reshape does not reserve new memory while EEG.data(:,:) does
end;	


com = sprintf( '%s = pop_iirfilt( %s, %s, %s, [%s], [%s]);', inputname(1), inputname(1), ...
			num2str( locutoff), num2str( hicutoff), num2str( trans_bw ), num2str( revfilt ) );
return;
