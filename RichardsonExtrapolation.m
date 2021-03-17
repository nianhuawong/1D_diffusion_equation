clc;clear;
U = [34.476349
34.476347
34.476346
];

u1 = U(1); 
u2 = U(2);
u3 = U(3);

h1 = 1/200; h2 = 1/100; h3 = 1/50;
e21 = u2 - u1;
e32 = u3 - u2;
r21 = h2 / h1;
r32 = h3 / h2;
s = sign(e32 / e21);

syms p;
q = log((r21^p-s)/(r32^p-s));

p=solve(p == 1.0/log(r21) * abs(log(e32 / e21)+q),p);
p=vpa(p,4)


