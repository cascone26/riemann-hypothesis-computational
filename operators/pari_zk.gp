default(realprecision, 20);
L = lfuncreate(x^2+1);
z = lfunzeros(L, 200);
print("COUNT ", #z);
for(i=1, #z, write("zk_zeros.txt", z[i]));
for(i=1, #z, my(d = lfun(L, 1/2 + I*z[i], 1)); write("zk_argLp.txt", arg(d)));
print("DONE");
