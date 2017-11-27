% Author :  by Serhat ÇAðÐDAÞ
% Reference to Entropy estimators - Improvements and Comparisons
% Wieczorkowski et al. (1999)
 


function h=vasicekm(v,m,n) 

    vals=sort(v);
    
    if(size(vals,1) == 1)
        vals = vals';
    end
  % Note that the intervals overlap for this estimator.
%   intvals = vals(2*m+1:n)-vals(1:n-2*m);
    intvals = [vals(2)-vals(1); vals(2*m+1:n)-vals(1:n-2*m); vals(n)-vals(n-1)];
%   n_int = n-2*m;
  
  hvec = log((n/(2*m))*intvals);
  h=(1/n)*sum(hvec);
