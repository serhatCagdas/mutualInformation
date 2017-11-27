function [ Y ,A] = generate_bivariate_uniform( mu, sigma, N )

    m = 1;
    
    X = (rand(N,2) - 0.5); 
    
    [V,D] = eig(sigma);
    
    A = (D.^0.5)*V*sqrt(12);
    Y=X*A;

%     Hj = 2*log(2) + log(abs(det((D.^0.5)/pinv(V)))) 
%     I = log(2)
end

