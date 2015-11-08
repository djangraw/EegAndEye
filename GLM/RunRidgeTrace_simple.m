function [betas,Xrss,Yrss] = RunRidgeTrace_simple(X,Y,lambdas)

% Created 2/25/15 by DJ.

p = size(X,2);
D = size(Y,2);
L = numel(lambdas);
betas = zeros(D,p,L);

% normalize Y
Yrss = nan(1,D);
for j=1:D
    Yrss(j) = sqrt(Y(:,j)'*Y(:,j));
    Y(:,j) = Y(:,j)/Yrss(j);
end
% normalize X
Xrss = nan(1,p);
for j=1:p
    Xrss(j) = sqrt(X(:,j)'*X(:,j));
    X(:,j) = X(:,j)/Xrss(j);
end


XTX = full(X'*X);    
XTY = full(X)'*Y;
% ridge trace
for i=1:numel(lambdas)
    betas(:,:,i) = ( (XTX+lambdas(i)*eye(p))^(-1)*XTY )'; % size Dxp            
end

% switch betas back to original units
betas = betas.*repmat(Yrss',1,p,L)./repmat(Xrss,D,1,L);