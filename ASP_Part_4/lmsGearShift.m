function [yhat, error, W] = lmsGearShift(x, z, mu, order, rho)
    N = length(x);
    yhat = zeros(1,N);
    error = zeros(1,N);
    W = zeros(order,N+1);
    MU = mu*ones(1,N);
    for i = order+1:N
        xpast = x(i-order+1:i);
        yhat(i) = W(:,i)'*xpast';
        error(i) = z(i)-yhat(i);
        W(:,i+1) = W(:,i) + MU(i)*error(i)*xpast';
        if error(i) <= error(i-1)
            MU(i+1) = MU(i)-rho*(error(i)-error(i-1));
        end        
    end
    W = W(:,2:end);
    
end