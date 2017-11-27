function [ H ] = uni_entr( a,b )
%transformed uniform joint entropy calculator 

 H = -(2/(a*b))*(a^2*log(a)/2 - a^2/4 - log(a*b)*a^2/2) +(b-a)/b*log(b); 


end

