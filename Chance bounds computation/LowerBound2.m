function LowerBound2

N = 250;
alpha = 0.15;
lossRate = 0.2;
Nalpha  = ceil(N*alpha);
rho = 0;
for i = 1 : Nalpha
    rho = rho + nchoosek(N, i)*lossRate^i*(1-lossRate)^(N-i);
end


display(rho)
end