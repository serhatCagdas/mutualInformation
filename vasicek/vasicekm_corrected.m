function [ hvec ] = vasicekm_corrected( v,m )

% Author      : Serhat ÇAÐDAÞ
% Description :Corrected version of mn spacing entropy estimator
% Reference   : Entropy estimators - Improvements and Comparisons
% Wieczorkowski et al. (1999)

n    = length(v);

Vmn  = vasicekm3(v,m,n);
n_int = n;% n-2*m;                  % n = m+1 ,......, n-m
hvec = Vmn - log(n_int) + log(2*m) - (1-2*m/n_int)*psi(2*m) + psi(n_int+1)...
       -(2/n_int)*sum(psi(1+m-1:2*m-1)); 
    
end

