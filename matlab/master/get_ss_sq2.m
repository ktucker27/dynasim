function sq2 = get_ss_sq2(n,o,f,g)

gamma = 1;

[~, szc, spc, smc] = coherent_matrix_onsite(n,o,f,g,gamma);
rhoss = coherent_rhoss(o, f, g, smc, spc);
[es, ess, ~] = calc_ev(rhoss, szc, spc, smc);
[sq2, ~, ~, ~] = get_squeezing(n, es, ess);

if abs(imag(sq2)) > 1e-6
    disp(['WARNING: Found large steady-state squeezing imaginary part: ', num2str(imag(sq2))]);
end

sq2 = real(sq2);