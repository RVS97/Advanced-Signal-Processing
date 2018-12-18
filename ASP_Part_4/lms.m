function [yhat, error, W] = lms(x, z, mu, order)
    N = length(x);
    yhat = zeros(1,N);
    error = zeros(1,N);
    W = zeros(order,N+1);
    for i = order+1:N
        xpast = x(i-order+1:i);
        yhat(i) = W(:,i)'*xpast';
        error(i) = z(i)-yhat(i);
        W(:,i+1) = W(:,i) + mu*error(i)*xpast';
    end
    W = W(:,2:end);
    
end