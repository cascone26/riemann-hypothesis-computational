default(realprecision, 25);
mf = mfinit([1,12,1], 0);
F = mfeigenbasis(mf)[1];
L = lfunmf(mf, F);
z = lfunzeros(L, 260);
print("COUNT ", #z);
for(i=1, #z, write("delta_zeros_pari.txt", z[i]));
print("DONE");
