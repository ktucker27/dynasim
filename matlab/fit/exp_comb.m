function val = exp_comb(b,t,t0)

t = t - t0;
val = zeros(size(t,1),1);
for i=1:2:size(b,1)
    val = val + b(i)*exp(1i*2*pi*b(i+1)*t);
end
