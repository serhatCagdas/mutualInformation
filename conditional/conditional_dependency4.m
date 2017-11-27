function [ I ] = conditional_dependency4( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
m = 1;
n_m = 1;
n_v = 1;
 
p = polyfit(X(:,1),X(:,2),n_m);
line_r_exp = polyval(p,X(:,1)); 

X2_p = (X(:,2) - line_r_exp); 
p = polyfit(X(:,1),X2_p.^2,n_v);
poly_var = polyval(p,X(:,1));
th = var(X2_p)*0.1;
poly_var = (poly_var >= th).*poly_var + (poly_var < th)*th;
poly_sd = sqrt(poly_var);

X2_p = X2_p./poly_sd;

H2  = vasicekm_corrected(X(:,2),m );
H2p = vasicekm_corrected(X2_p  ,m );
% 
% Hjoint = H1 + H2p + mean(log(poly_sd));
I  = H2 - H2p - mean(log(poly_sd));
end

