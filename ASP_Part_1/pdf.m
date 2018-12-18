function [values, centers] = mypdf(x)
    [counts, centers] = hist(x,100);
    values = counts/trapz(centers,counts);
end