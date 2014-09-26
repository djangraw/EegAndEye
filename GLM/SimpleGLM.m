function [h,tH,Efit,Eres,lambda] = SimpleGLM(E,t,S,regressor_range,reg_method,lambda,demeanX)

% Implements a simple General Linear Model to EEG data as suggested in 
% Dandekar et al., J. Neurosci., 2011.
%
% [h,tH,Efit,Eres] = SimpleGLM(E,t,S,regressor_range,reg_method,lambda,demean)
%
% INPUTS:
% - E is an LxD matrix of EEG data, where L is the # of samples and D is 
%  the # of channels.
% - t is a column vector of the corresponding times, also of length L.
% - S is an L x Nr*(Nh+1) matrix, where Nr is the number of regressors
%  and Nh is the number of samples in each desired response function.  It 
%  specifies the regressor event times at various offsets and is created 
%  by the program GetGlmRegressors.m. 
% - regressor_range is a scalar or 2-element vector indicating the number 
%  of samples on either side of each regressor event that are assumed to be
%  affected by a regressor event.  This should be the same number used as 
%  an input to GetGlmRegressors when making S. 
% - reg_method is a string indicating the type of regression to use.  The
%  options are:
%   + 'mvregress', simple multivariate linear regression (default)
%   + 'diagonal', like mvregress but the covariance matrix is forced to be
%   diagonal.
%   + 'glmfit', matlab's generalized linear model function.  (HASN'T BEEN
%   SUCCESSFULLY USED YET!  NEEDS WORK!)
%   + 'leastsquares', simple Maximum Likelihood solution. ONLY THIS OPTION
%   SUPPORTS MULTIPLE CHANNELS!
%   + 'ridge', ridge regression in which a diagonal matrix lambda*I is
%   added to the covariance matrix to help deal with non-orthogonal
%   regressors. Needs lambda input. SUPPORTS MULTIPLE CHANNELS! 
%   [default = 1e-6]. 
% - demeanX is a scalar indicating whether you want to de-mean the regressor
%   matrix (in the case of ridge regression).
%
% OUTPUTS:
% - h is an Nr x Nh x D matrix in which each row i, page j is the inferred 
%  response of the EEG to event i on channel j.
% - tH is a vector of length Nh indicating the corresponding times
%  relative to the event.
% - Efit is a vector of length L containing the EEG data as reconstructed
%  by the algorithm.
% - Eres is a vector of length L containing the residuals (E-Efit).
%
% Created 12/22/11 by DJ.
% Updated 1/10/12 by DJ - added other regression options (use ridge if your
% regressors are not orthogonal)
% Updated 1/12/12 by DJ - comments
% Updated 3/22/13 by DJ - changed input from Nt to regressor_range, allow
%  2-element vector.
% Updated 4/4/13 by DJ - fixed Nt holdover bug
% Updated 8/3/14 by DJ - Added least squares option, multi-channel support
% Updated 8/8/14 by DJ - added ridge option, lambda input
% Updated 8/19/14 by DJ - standardized S columns for ridge
% Updated 8/26/14 by DJ - made de-meaning optional in ridge

if nargin<1 || isempty(E)
    E = [0 0 0 1 -2 1 0 -2 0]'; % column vector
end
if nargin<2 || isempty(t)
    t = 1:size(E,1);
end
if nargin<3 || isempty(S)
    [~,S] = GetGlmRegressors_v2p0(t,{5,[5 8]},[],regressor_range); % test case
end
if nargin<4 || isempty(regressor_range)
    regressor_range = 0; % only use exact event time
end
if nargin<5 || isempty(reg_method)
    reg_method = 'mvregress';
end
if nargin<6 || isempty(lambda)
    lambda = NaN;
end
if nargin<7 || isempty(demeanX)
    demeanX = false;
end

% Handle 1-element range
if numel(regressor_range)==1
    regressor_range = [-regressor_range regressor_range];
end


% Set options
plotResults = 0;

% Get constants
% L = size(E,1); % Length of EEG data
D = size(E,2);
Nh = diff(regressor_range)+1; % length of response vectors
% Nh = 2*Nt+1; % length of response vectors
dt = median(diff(t)); % time between samples
Nr = size(S,2)/Nh; % # regressors

% Get full S matrix
S = full(S);

% Response matrix
switch reg_method
    case 'spm'
        fprintf('spm...')
        results = spm_glm(E,double(S));
        H = results.w;
        % Find fit and residuals
        if nargout>2 || plotResults
            Efit = S*H;
            Eres = E-Efit;
        end
    case 'spm_robust'
        fprintf('spm_robust...')
        [~,H] = spm_robust_glm(E,double(S),1,3);
        % Find fit and residuals
        if nargout>2 || plotResults
            Efit = S*H;
            Eres = E-Efit;
        end
    case 'glmfit'
        fprintf('glmfit...')
        H = glmfit(S,E,'normal','constant','off');
        % Find fit and residuals
        if nargout>2 || plotResults
            Efit = S*H;
            Eres = E-Efit;
        end
%     case 'ridge'
%         fprintf('ridge regression...')
%         H = ridge(E,S,1e-6); % 3rd input is small regularization term
%         % Find fit and residuals
%         if nargout>2 || plotResults
%             Efit = S*H;
%             Eres = E-Efit;
%         end
    case 'mvregress'
        fprintf('using mvregress...')
        [H,~,Eres] = mvregress(S,E);
        % H = pinv(S'*S)*S'*E; % Equation 3
        % Find fit and residuals
        Efit = E-Eres;
    case 'diagonal' % mvregress with diagonal covariance matrix
        fprintf('diagonal cov matrix...')
        [H,~,Eres] = mvregress(S,E,'covtype','diagonal');
        % H = pinv(S'*S)*S'*E; % Equation 3
        % Find fit and residuals
        Efit = E-Eres;
    case 'leastsquares'
        fprintf('least squares solution...')
        H = (S'*S)^(-1)*S'*E; % try it!
        Efit = S*H;
        Eres = E-Efit;
    case 'ridge'                
        if demeanX
            fprintf('De-meaning and ');
            mnS = mean(S,1);
        else
            fprintf('NOT de-meaning, but ');
            mnS = zeros(1,size(S,2));
        end
        fprintf('Standardizing S matrix...\n')
        parfor i=1:size(S,2)            
            S(:,i) = (S(:,i) - mnS(i)) / sqrt(sum((S(:,i) - mnS(i)).^2));
        end
        if isequal(lambda,'auto')
            fprintf('starting with least squares solution...\n')
            H0 = (S'*S)^(-1)*S'*E; % try it!            
            fprintf('using results to calculate lambda...\n')
            Eres0 = E-S*H0;
            [n,p] = size(S);
            lambda = nan(1,D);
            fprintf('Ridge Regression:\n');
            for j=1:D
                sigmasq = (Eres0(:,j)'*Eres0(:,j))/(n-p);
                lambda(j) = p*sigmasq/sum(H0(:,j).^2);            
            end
            H = nan(p,D);
            STS = S'*S;
            STE = S'*E;
            parfor j=1:D
                tic
                fprintf('Electrode %d/%d, lambda=%.3g...',j,D,lambda(j))
%                 H(:,j) = (S'*S + lambda(j)*eye(size(S,2)))^(-1)*S'*E(:,j); % single elec, single lambda.
                H(:,j) = (STS + lambda(j)*eye(p))^(-1)*STE(:,j); % single elec, single lambda.
                t = toc;
                fprintf('Done! took %.2f seconds.\n',t);               
            end
        else
            fprintf('ridge regression, lambda=%.3g...',lambda)
            H = (S'*S + lambda*eye(size(S,2)))^(-1)*S'*E; % All elecs with same lambda
        end
        Efit = S*H;
        Eres = E-Efit;
    otherwise
        error('regression method %s not recognized!',reg_method);
end

% Response functions
h = zeros(Nr,Nh,D);
for j=1:D
    for i=1:Nr
        h(i,:,j) = H((i-1)*Nh+(1:Nh),j);
    end
end
tH = (regressor_range(1):regressor_range(2))*dt; % corresponding time vector, assuming constant sampling rate

% Plot results
if plotResults
    if D>1
        disp('Plotting >1 electrode not yet supported.')
        return;
    end
    figure(223)
    plot(tH,h') % plot responses
    xlabel('time (same units as input t)')
    ylabel('EEG response (same units as input E)')
    title('EEG response to regressors')
    names = cell(1,Nr);
    for i=1:Nr
        names{i} = sprintf('Regressor %d',i);
    end
    legend(names);

    % Get regressors
    s = S(:,1-regressor_range(1):Nh:end)'; % events at t=0 points
    figure(224)
    handles = nan(1,Nr+1);
    handles(1) = subplot(Nr+1,1,1);
    plot(t,[E,Efit]);
    legend('true EEG','Reconstructed EEG');
    ylabel('EEG (same units as input E)');
    for i=1:Nr
        handles(i+1) = subplot(Nr+1,1,i+1);
        plot(t,s(i,:));
        ylabel(sprintf('Regressor %d',i))    
    end
    xlabel('time (same units as input t)')
    linkaxes(handles,'x');
end