function CHuge = MakeRepeatedImage(levelname,offset,scale,repeatoffset,nRepeats)

% Creates a repeated version of an image by adding it to itself with an
% offset.
%
% CHuge = MakeRepeatedImage(levelname,offset,scale,repeatoffset,nRepeats)
%
% INPUTS:
% -levelname is a string that indicates an image filename in the current
% path.
% -repeatoffset is a 2-D vector indicating the x and y offset for each 
% repeat of the image.
% -nRepeats is a scalar indicating the number of desired repeats.
%
% OUTPUTS:
% -CHuge is an image matrix containing the repeated image.  It can be saved 
% using imwrite.
%
% Created 4/15/11 by DJ.

% load image
C = imread(levelname);
% plot image
subplot(1,2,1)
imagesc(C);
axis image
hold on;
% overlay scale bars
plot([offset(1),offset(1)], [offset(2), offset(2)+scale(2)*100],'g');
plot([offset(1),offset(1)+scale(1)*100], [offset(2), offset(2)],'g');

% make big canvas
CHuge = uint8(zeros(size(C,1)+round(abs(repeatoffset(1)*scale(1)*nRepeats)+1), size(C,2)+round(abs(repeatoffset(2)*scale(2)*nRepeats)+1),3));
xStart = size(CHuge,1) - size(C,1);
yStart = size(CHuge,2) - size(C,2);
% add repeats
for i=0:nRepeats
    xToChange = xStart - round(repeatoffset(1)*scale(1)*i)+(1:size(C,1));
    yToChange = yStart - round(repeatoffset(2)*scale(2)*i)+(1:size(C,2));
    
    CHuge(xToChange,yToChange,:) = CHuge(xToChange,yToChange,:) + C;
end
% plot repeated image
subplot(1,2,2)
imagesc(CHuge);
axis image
hold on;
% overlay scale bars
plot([offset(1),offset(1)], [offset(2), offset(2)+scale(2)*100],'g');
plot([offset(1),offset(1)+scale(1)*100], [offset(2), offset(2)],'g');

