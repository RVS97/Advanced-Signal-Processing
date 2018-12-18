function [xhat, error, A] = lmsSignRegressor(x, mu, order, rho)
    N = length(x);
    xhat = zeros(1,N);
    A = zeros(order, N+1);
    error = zeros(1,N);
    MU = mu*ones(1,N);
    for i = order+1:N
        xpast = x(i-1:-1:i-order);
        xhat(i) = A(:,i)'*xpast';
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + MU(i)*error(i)*sign(xpast');
        if error(i) <= error(i-1)
            MU(i+1) = MU(i)-rho*(error(i)-error(i-1));
        end        
    end
    A = A(:,2:end);
    
end