default(realprecision, 20);
mf = mfinit([23, 1, -23], 1);
F = mfeigenbasis(mf)[1];
L = lfunmf(mf, F);
z = lfunzeros(L, 200);
print("COUNT ", #z);
for(i=1, #z, write("h23_zeros.txt", z[i]));
for(i=1, #z, my(d = lfun(L, 1/2 + I*z[i], 1)); write("h23_argLp.txt", arg(d)));
print("DONE");
