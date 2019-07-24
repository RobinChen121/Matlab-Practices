function [ kurt] = Kurt(d, p)

mean = Mean(d, p);
variance = Variance(d, p);
kurt = (d - mean).^4 * p'/ variance^2;


end

