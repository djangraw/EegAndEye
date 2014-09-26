function [data_sorted, order] = ImageSortedData(data,x,y,rt,mode)

% Images the given data with rows re-ordered according to a sorted vector.
%
% [data_sorted, order] = ImageSortedData(data,x,y,rt,mode)
%
% INPUT:
% -data is an MxN matrix of the data you want to image.
% -x is an N-element vector with the x axis of your plot
% -y is an M-element vector with the y axis of your plot
% -rt is an M-element vector of the values you want to use to sort the data
%  (e.g., reaction times)
% -mode is a string used as an input to maltb's sort function 
%  ('ascend' (default) or 'descend')
%
% OUTPUT:
% -data_sorted is the data matrix with rows re-ordered by rt
% -order is the order of the rows (i.e., data_sorted = data(order,:);
%
% Created 9/4/12 by DJ.
% Updated 9/12/12 by DJ - added y input, outputs.
% Updated 10/8/12 by DJ - added mode input

if nargin<5
    mode = 'ascend';
end

% Sort
[rt_sorted, order] = sort(rt,mode);
data_sorted = data(order,:);

% Plot
hold on;
imagesc(x,y,data_sorted);
plot([0 0],[min(y) max(y)],'k-')
plot(rt_sorted,y,'k-','linewidth',2);
