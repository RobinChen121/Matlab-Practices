function LowerBound

N = 1250;
r = 260;
p  = 0.95
tempP = 0
lossRate = 0.228
%while tempP <= 1 - p
for i = 1 : r
    tempP = tempP + nchoosek(N, i)*lossRate^i*(1-lossRate)^(N-i);
end
%end

display(tempP)
end