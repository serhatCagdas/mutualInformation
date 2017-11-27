function [ I ] = kraskov_MI_light( X )
k = 1;

N = size(X,1);
nx       = zeros(N,1);
ny       = zeros(N,1);
epsx     = zeros(N,1);
epsy     = zeros(N,1); 
eps      = zeros(N,1); 

for i =  1: N  
         
     d = abs(X - repmat(X(i,:),N,1)); 
     dist_vec = max(d,[],2);%sqrt(sum((d.*d),2));
     [dist_v_sotrted, ind] = sort(dist_vec, 'ascend');
     epsx(i) = max(d(ind(1:k+1),1));
     epsy(i) = max(d(ind(1:k+1),2));
     eps(i)  = dist_v_sotrted(k+1); 

     nx(i) = sum(d(:,1) < eps(i)) - 1;
     ny(i) = sum(d(:,2) < eps(i)) - 1;
               
end
 
    I = psi(k) - mean(psi(nx+1)+ psi(ny+1)) + psi(N);
  
 
end

