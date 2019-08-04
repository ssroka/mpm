function int = integrate_logspace(X,Y)
% Numerical integration over a logarithmic space.
% The function approximates the integral with samples X
% that are spaced logarithmically.
% The integration is accurate where the function Y is piece-wise
% linear over a logaritchmic scale.
% When the function is smooth on a logarithmic scale,
% the integral may be evaluated accurately, using fewer
% samples in comparison to linear sampling, and the
% integral is evaluated much faster.
% A natural application is integration of transfer functions,
% for which the frequency samples are usually logarithmically
% spaced.
%
% Inputs:
% X - vector of samples.
% Y - corresponding function values. If Y is a matrix,
% each of its columns is integrated with respect to X.
% Output:
% int - integration result. If Y is a matrix, int is a row
% vector which elements correspond to integrals of Y columns.
%
% Handling zero values in Y:
% a) if a column of Y is all zeros the corresponding output is zero.
% b) if a column of Y contains both zero and non-zero elements the output
% is NaN
%
% Written by Dr. Yoash Levron, Dec. 2015

int = NaN;
NN = length(X); % number of samples
X=reshape(X,NN,1); % X is now a column vector
if isvector(Y)
    Y=reshape(Y,length(Y),1); % Y is now a column vector
end
[tt, R] = size(Y); % R - number of functions to integrate
if (tt~=NN)
    disp('Error in ''integrate_logspace'' - input/output dimensions must agree');
    return
end

col_zeros = zeros(1,R);  % indexes of columns of Y containing zeros
for kk = 1:R
    col_zeros(kk) = isempty(find(Y(:,kk)));
end
col_zeros_ind = find(col_zeros); % indices of columns which are all zeros
col_nonzeros_ind = find(~col_zeros); % indices of columns which are not all zeros

% rearrange input to contain only non-zero columns:
orgR = R;
if (~isempty(col_zeros_ind))
    Y = Y(:,col_nonzeros_ind);
    R = size(Y,2);
end

% compute integral
eps0 = 1e-15;
if (R==1)
    logX = log(X);
    logY = log(Y);
    m = diff(logY)./diff(logX);
    n = logY(1:(end-1)) - m.*logX(1:(end-1));
    mp1 = m+1;
    p = (exp(n)./mp1).*(X(2:end).^mp1 - X(1:(end-1)).^mp1);
    ind = find(abs(mp1)<eps0);
    p(ind) = exp(n(ind)).*(log(X(ind+1)./X(ind)));
    int_a = sum(p);
else   % R>1
    logX = log(X);
    logY = log(Y);    
    m = diff(logY)./repmat(diff(logX),1,R);
    n = logY(1:(end-1),:) - m.*repmat(logX(1:(end-1)),1,R);
    mp1 = m+1;
    X2 = repmat(X(2:end),1,R);
    X1 = repmat(X(1:(end-1)),1,R);
    p = (exp(n)./mp1).*(X2.^mp1 - X1.^mp1);
    ind = find(abs(mp1)<eps0);
    [I, J] = ind2sub(size(m),ind);
    p(ind) = exp(n(ind)).*(log(X(I+1)./X(I)));
    int_a = sum(p);    
end

% handle all-zero columns:
int = zeros(1,orgR);
int(col_nonzeros_ind) = int_a;

end



