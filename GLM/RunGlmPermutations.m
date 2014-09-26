function [b_data,b_rand,iPerm] = RunGlmPermutations(XXinvX,Y,nPerm,Nh)

% Calculate
disp('Calculating response functions...')
B_data = XXinvX*Y;

% Arrange into matrix
Nr = size(XXinvX,1)/Nh;
D = size(Y,2);
b_data = zeros(Nr,Nh,D);
for j=1:D
    for k=1:Nr
        b_data(k,:,j) = B_data((k-1)*Nh+(1:Nh),j);
    end
end
b_data = permute(b_data,[3 2 1]); % elec x time x regressors
%% Yfit = X*B;
% Yres = Y-Yfit;

[b_rand,iPerm] = deal(cell(1,nPerm));
for i=1:nPerm
    % Calculate
    fprintf('Permutation %d/%d...\n',i,nPerm)
    iPerm{i} = randperm(size(Y,1));
    B_rand = XXinvX*Y(iPerm{i},:); % try it!  
    % Arrange into matrix
    b_rand{i} = zeros(Nr,Nh,D);
    for j=1:D
        for k=1:Nr
            b_rand{i}(k,:,j) = B_rand((k-1)*Nh+(1:Nh),j);
        end
    end
    b_rand{i} = permute(b_rand{i},[3 2 1]); % elec x time x regressors
    
end
disp('Done!')