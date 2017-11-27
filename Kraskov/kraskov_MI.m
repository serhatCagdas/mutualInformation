function [ I ] = kraskov_MI( X , k, is_common_eps)
%UNTÝTLED4 Summary of this function goes here
%   Detailed explanation goes here
alpha = 0.2;
if k == 3
    alpha = 0.13;
elseif(k == 4)
    alpha = 0.2;
elseif k==5
    alpha = 0.3;
end

N = size(X,1);


nx       = zeros(N,1);
ny       = zeros(N,1);
epsx     = zeros(N,1);
epsy     = zeros(N,1); 
eps      = zeros(N,1); 

for i =  1: N  
         
             d = abs(X - repmat(X(i,:),N,1)); 
             dist_vec = max(d');%sqrt(sum((d.*d),2));
             [dist_v_sotrted, ind] = sort(dist_vec, 'ascend');
             epsx(i) = max(d(ind(1:k+1),1));
             epsy(i) = max(d(ind(1:k+1),2));
             eps(i)  = dist_v_sotrted(k+1);
             
             if is_common_eps == 1
                 nx(i) = sum(d(:,1) < eps(i)) - 1;
                 ny(i) = sum(d(:,2) < eps(i)) - 1;
             elseif is_common_eps == 3
                 nx(i) = sum(d(:,1) <= epsx(i)) -1;
                 ny(i) = sum(d(:,2) <= epsy(i)) -1; 
                 coeff = pca(X(ind(1:k+1),:));
                 lc1 = X(ind(1:k+1),:)*coeff(:,1);
                 lc2 = X(ind(1:k+1),:)*coeff(:,2);
                 Vm  = (max(lc1)-min(lc1))*(max(lc2)-min(lc2));
                 V = (2*epsx(i))*(2*epsy(i));
                 if (Vm/V < alpha)
                    LNC(i) = log10(Vm/V);
                 else
                    LNC(i) = 0;
                 end
                 
             else
                 nx(i) = sum(d(:,1) <= epsx(i)) -1;
                 ny(i) = sum(d(:,2) <= epsy(i)) -1; 
             end 
             
             
end


 if is_common_eps == 1
    I = psi(k) - mean(psi(nx+1)+ psi(ny+1)) + psi(N);
 elseif is_common_eps == 3
    I = psi(k) - 1/k - mean(psi(nx)+ psi(ny)) + psi(N) - mean(LNC);
 else
    I = psi(k) - 1/k - mean(psi(nx)+ psi(ny)) + psi(N);
 end 

 
%  if I < 0
%      I = 0;
%  end 
end

