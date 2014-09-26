% MakeGrating.m
%
% Plots a high-contrast and a low-contrast sinusoid at the angle and period
% specified.  These can be used as textures for the stimuli in Unity.
%
% Created 4/5/11 by DJ.

% set params
gratesize = 200; % size of image in pixels
period = 50; % period of grating in pixels
angle = 90; % angle relative to horizontal in degrees
highcontrast = 1;
lowcontrast = .5;

% create building blocks
A = zeros(gratesize);
[X,Y] = meshgrid(1:gratesize);
Z = X*cos(angle*pi/180)+Y*sin(angle*pi/180);

% create high-contrast grating
B = highcontrast*sin(2*pi*Z/period);
subplot(1,2,1)
imagesc(B)
set(gca,'CLim',[-1 1],'xtick',[],'ytick',[]);

% create low-contrast grating
C = lowcontrast*sin(2*pi*Z/period);
subplot(1,2,2)
imagesc(C)
set(gca,'CLim',[-1 1],'xtick',[],'ytick',[]);