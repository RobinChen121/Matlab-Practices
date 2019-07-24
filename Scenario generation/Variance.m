function [ variance ] = Variance( d, p )

mean = Mean(d, p);

T = length(d);
variance = 0;
for i = 1: T
    variance = variance + p(i)*(d(i) - mean)^2;
end
% variance = (d - mean).^2 * p';

end

