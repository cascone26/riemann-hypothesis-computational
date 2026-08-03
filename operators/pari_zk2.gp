default(realprecision, 25);
L = lfuncreate(x^2+1);
z = lfunzeros(L, 200, 32);
print("COUNT ", #z);
for(i=1, #z, write("zk_zeros_v2.txt", z[i]));
print("DONE");
