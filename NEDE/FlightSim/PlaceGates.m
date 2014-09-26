function PlaceGates(d,N,gate_type,params)
%
% INPUTS:
% -d is a scalar indicating the distance between gates, in meters
% -N is a scalar indicating number of gates
% -gate_type is a string ('cosine' and 'steps' are currently supported)
% -params is a struct containing info about the path, including frequencies
% and amplitudes, widths, etc.
%
% Created 9/3/14 by DJ.

%% Declare gate parameters
if ~exist('d','var') || isempty(d)
    d = 500; % dist btw gates, in meters
end
if ~exist('N','var') || isempty(N)
    N = 30; % # of gates
end
if ~exist('gate_type','var') || isempty(gate_type)
    gate_type = 'steps';
end
% set positions of gates
zpos = 0:d:(d*(N-1)); % in meters


if nargin<4
    auto_params = true;
else
    auto_params = false;
end


switch gate_type
    case 'cosine'
        %% Cosine wave gates

        if auto_params
            x_period = 10000; % in meters
            x_amp = 0;
            x_offset = 0;
            y_period = 8000; % in meters
            y_amp = 100;
            y_offset = 150;
        else
            UnpackStruct(params);
        end
        % set positions
        xpos = x_offset - x_amp/2*cos(2*pi*zpos/x_period);
        ypos = y_offset - y_amp/2*cos(2*pi*zpos/y_period);
        % set widths
        width = repmat(50,size(zpos));
        filename = sprintf('CosinePath_amp%d-%d_pd%d-%d_off%d-%d.txt',x_amp,y_amp,x_period,y_period,x_offset,y_offset);
    case 'steps'
        if auto_params
            height = 100;
            widthRange = [150 0];
            widthStep = 30;
        else
            UnpackStruct(params)
        end
        % set positions
        xpos = zeros(size(zpos));
        ypos = repmat(height,size(zpos));
        % set widths
        width = ceil(linspace(widthRange(1),widthRange(2),N)/widthStep)*widthStep;
        filename = sprintf('Steps_range%d-%d_step%d.txt',widthRange,widthStep);
    otherwise
        error('gate type not recognized!');        
        
end

%% Write
% write results to file
fileID = fopen(filename,'w');
for i=1:N
    fprintf(fileID,'%.2f, %.2f, %.2f, %.2f\n',xpos(i),ypos(i),zpos(i),width(i));
end
fclose(fileID);

%% Plot
clf
subplot(3,1,1);
hold on
% plot path
plot(zpos,xpos,'r.-');
% Plot rings
for i=1:numel(zpos)
    plot([zpos(i), zpos(i)], [xpos(i)-width(i)/2, xpos(i)+width(i)/2],'r-');        
end
plot(zpos,xpos-width/2,'r:');
plot(zpos,xpos+width/2,'r:');
xlabel('z')
ylabel('x');
title('view from above')

subplot(3,1,2);
hold on
% plot path
plot(zpos,ypos,'r.-');
% Plot rings
for i=1:numel(zpos)
    plot([zpos(i), zpos(i)], [ypos(i)-width(i)/2, ypos(i)+width(i)/2],'r-');        
end
plot(zpos,ypos-width/2,'r:');
plot(zpos,ypos+width/2,'r:');
xlabel('z')
ylabel('y');
title('view from the side')

subplot(3,1,3);
plot3(zpos,xpos,ypos,'r.-');
set(gca,'View',[-5.5 56]);
xlabel('z')
ylabel('x');
zlabel('y');
title('3D view')

MakeFigureTitle(filename);


