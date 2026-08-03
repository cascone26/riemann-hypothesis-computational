default(realprecision, 20);
mf = mfinit([1,12,1], 0);
F = mfeigenbasis(mf)[1];
L = lfunmf(mf, F);
z = lfunzeros(L, 600);
print("COUNT ", #z);
for(i=1, #z, write("delta_zeros_600.txt", z[i]));
for(i=1, #z, my(d = lfun(L, 6 + I*z[i], 1)); write("delta_argLp_600.txt", arg(d)));
print("DONE");
