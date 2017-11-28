function [ I ] = conditional_dependency( X,opt)
%AUTHOR - Serhat CAGDAS - IZTECH 2017
%@param 
% X : Nx2 data to find mutual information
% opt : 1-> polynomal fit
%       2-> convolutional fit

MAX_DEG  = 10;       % max polynomial degree
MAX_LN   = 9;       % max cov len 
LEN_STEP = 0.2;     % conv len step size
 
N = length(X(:,1));
m = 1;

H1  = vasicekm_corrected(X(:,1),m );
H2  = vasicekm_corrected(X(:,2),m );
 

if (opt == 1)       % polynomial fit       
    Hj   = zeros(MAX_DEG,MAX_DEG);  
    
    for d_m = 1: MAX_DEG

         %X2_p = x2 - E[X2| X1 = x1]
         p = polyfit(X(:,1),X(:,2),d_m);  % mean fit
         mean_fit = polyval(p,X(:,1));  
         X2_p = (X(:,2) - mean_fit);

        for d_v = 1:MAX_DEG

            % X2_p = sqrt(x2p /(var(X2_p | X1 = x1)) )
            p = polyfit(X(:,1),X2_p.^2,d_v); % variance fit
            var_fit = polyval(p,X(:,1));

            th = var(X2_p)*0.1;                         % variance can't be negative
            var_fit = (var_fit >= th).*var_fit + (var_fit < th)*th;
            sd_fit = sqrt(var_fit); 
            X2_pp = X2_p./sd_fit;

            % joint entropy is H(X1) + H(X2') + mean(log(sd)) 
            H2p = vasicekm_corrected(X2_pp  ,m ); 
            Hj(d_m,d_v) = H1 + H2p + mean(log(sd_fit));

        end 
    end 

else                %convolutional fit
    Hj = zeros(MAX_LN,MAX_LN); 
    
    [X1,ind] = sort(X(:,1),'ascend');   %sort data
     X2 = X(ind,2);
     X = [X1 X2];

    mean_arr = fibonacci(2:10)/fibonacci(11);
    
    for d_m = 1: MAX_LN
        mcl = ceil(N*mean_arr(d_m)*LEN_STEP);
        mcl = mcl + (mcl == 1); 
        
        %X2' = X2 - E[X2| X1 = x1]
        num_of_el = conv(ones(N,1),ones(mcl,1),'same'); % mean fit
        mean_fit  = conv(X(:,2) , ones(mcl,1),'same')./num_of_el;
        X2_p = (X(:,2) - mean_fit);

        for d_v = 1: MAX_LN
            vcl = ceil(N*mean_arr(d_v)*LEN_STEP);
            vcl = vcl + (vcl == 1);
            % X2' = sqrt(x2p /(var(X2_p | X1 = x1)) )
            num_of_el2 = conv(ones(N,1),ones(vcl,1),'same');% variance fit
            var_fit   = conv(X2_p.^2 , ones(vcl,1),'same')./num_of_el2;
            sd_fit = sqrt(var_fit); %std dev
            X2_pp = X2_p./sd_fit;
            % joint entropy is H(X1) + H(X2') + mean(log(sd)) 
            H2p = vasicekm_corrected(X2_pp  ,m ); 
            Hj(d_m,d_v) = H1 + H2p + mean(log(sd_fit));

        end
    end 
end             % end of opt if
 
I =H1 + H2 - min(min(Hj));
end