function skew = Skew(d, p)

mean = Mean(d, p);
variance = Variance(d, p);
skew = (d - mean).^3 * p'/ variance^1.5;
end