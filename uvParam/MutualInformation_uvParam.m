function I = MutualInformation_uvParam(X)

x = X(:,1);
y = X(:,2);
m = 1;
N = 181; % divide the [0 pi] into N equal-sized bins...

% Dimension Correction...
if size(x, 1) > size(x, 2)
    X = x';
else
    X = x;
end

if size(y, 1) > size(y, 2)
    Y = y';
else
    Y = y;
end

% Initialization...
P = [X; Y];
theta = linspace(-90, 89, N);
U = @(phi) [cosd(phi); sind(phi)]; u = U(theta);
Ht = zeros(N, 1);

% Marginal Entropy Estimation...
for n = 1:N
    % Enter Your Own Entropy Code...
    Ht(n) = vasicekm_corrected(u(:, n)' * P,m);
    %Ht(n) = KL_Marginal(u(:, n)' * P);    
end

% Matrix Operations for calculation of Joint Entropy...
Hy = tril(ones(N-1) .* Ht(2:end));
Hx = triu(ones(N-1) .* Ht(1:end-1))';
Hk = (Hy + Hx) + (Hy + Hx)' - diag(diag(Hy + Hx));

Dy = tril(ones(N-1) .* theta(2:end)');
Dx = triu(ones(N-1) .* theta(1:end-1)')';
Dk = (Dy - Dx) + (Dy - Dx)' - diag(diag(Dy - Dx));

Hij = Hk - log(abs(sind(Dk)));
[v,i] = min(Hij);
[v,j] = min(v);
I = vasicekm_corrected(x,m) + vasicekm_corrected(y,m) - min(Hij(:));

% End Of Routine...
end