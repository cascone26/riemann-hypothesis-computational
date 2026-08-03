default(realprecision, 20);
mf = mfinit([1,12,1], 0);
F = mfeigenbasis(mf)[1];
L = lfunmf(mf, F);
z = readvec("delta_zeros_gp_input.txt")[1];
for(i=1, #z, my(d = lfun(L, 6 + I*z[i], 1)); write("delta_argLp_pari.txt", arg(d)));
print("DONE ", #z);
