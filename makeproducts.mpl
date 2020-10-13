#	We want to use this to generate the extra products 

restart;
with(Physics):
read `sta.mpl`:
with(sta);
read `lasenby-ashdown/dirac_maple_general_h_functions_v6.mpl`:
read `ewgt.mpl`:
with(ewgt);
read `tools.mpl`:
with(tools);


a:=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16];

b:=[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16];
writeto("sta_products_save.mpl");
lprint(ds(a&@b));
lprint(ds(a&.b));
lprint(ds(a&^b));
lprint(ds(a&!b));
writeto(terminal);
print(`done`);
