with(plots):
with(FileTools):
with(PDEtools):
with(DEtools):
with(StringTools):
with(CLIo):
read `tools.mpl`:
with(tools):
read `theory_tools.mpl`:

read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:

dit({},"First we would like to check (88) for Class 3C:");

for ii in H,q,X,Y do
  tmp[At][ii]:=simplify(subs(A=0,S3=0,S2=-1/3,S1=-1/3,U2=-4/3,Q(t)=w(t)/sqrt(K),tmp[Friedmann_form][ii])):
  #tmp[At][ii]:=simplify(subs(A=0,S3=0,S2=-1/3,S1=-1/3,U2=-4/3,U1=-K*l/3,Q(t)=w(t)/sqrt(K),tmp[Friedmann_form][ii])):
end do:

K:=1:

for ii in H,q,Y do
  tmp[aa][ii]:=simplify(subs(isolate(tmp[At][X],U(t)),tmp[At][ii])):
end do:

eprint(tmp,aa,"General definition of Class A");


substitutions:=w(t)=w0,H(t)=H,O_d(t)=0,O_r(t)=0,O_L(t)=L/(3*H^2):

for ii in H,q,Y do
  tmp[late][ii]:=simplify(subs(substitutions,tmp[aa][ii])):
end do:

eprint(tmp,late,"here is the frozen form");

sys:={tmp[late][H],tmp[late][q],tmp[late][Y]}:

solve(sys,{H,w0});

for ii in H,q,Y do
  tmp[ke][ii]:=numer(lhs(tmp[late][ii])-rhs(tmp[late][ii])):
end do:

sys:={tmp[ke][H],tmp[ke][q],tmp[ke][Y]}:

solve(sys,{H,w0});

print(%);

substitutions:=w(t)=w0,H(t)=Hz,O_d(t)=0,O_r(t)=0,O_L(t)=L/(3*Hz^2),U(t)=3*(f0+Hz):

for ii in H,q,Y,X do
  tmp[later][ii]:=simplify(subs(substitutions,tmp[At][ii])):
end do:

eprint(tmp,later,"here is the frozen form");

sys:={tmp[later][H],tmp[later][q],tmp[later][Y],tmp[later][X]}:

solve(sys,{H,w0,f0});

for ii in H,q,Y,X do
  tmp[ker][ii]:=numer(lhs(tmp[later][ii])-rhs(tmp[later][ii])):
end do:

sys:={tmp[ker][H],tmp[ker][q],tmp[ker][Y],tmp[ker][X]}:

solve(sys,{Hz,w0,f0});

print(%);

fin();

for ii in H,q,Y do
  tmp[b][ii]:=simplify(subs(substitutions,tmp[aa][ii])):
end do:

eprint(tmp,b,"substituting for matter");

for ii in H,q,Y do
  tmp[c][ii]:=simplify(subs(H(t)=diff(a(t),t)/a(t),tmp[b][ii])):
end do:

eprint(tmp,c,"substituting for Hubble");

for ii in H,q,Y do
  tmp[d][ii]:=simplify(dchange(t=tfun(x),tmp[c][ii],{x},known={tfun})):
  tmp[d][ii]:=evala(simplify(subs(diff(tfun(x),x)=a(x)/H_0,tmp[d][ii]))):
end do:

eprint(tmp,d,"switching to conformal time");

for ii in H,q,Y do
  tmp[e][ii]:=simplify(numer(normal(lhs(tmp[d][ii])-rhs(tmp[d][ii])))):
  tmp[e][ii]:=simplify(subs(O_l0=l,O_L0=L,O_m0=M^2,O_r0=R^2,tmp[e][ii])):
end do:

H_0:=1:

eprint(tmp,e,"tidying up");

for ii in H,q,Y do
  tmp[f][ii]:=primpart(tmp[e][ii]):
end do:

eprint(tmp,f,"more tidying up");

#(*
dit({},"power series for Class A");
w_s:=0:
a_s:=0:
for ii from 0 to 10 do
  w_s:=w_s+w||ii*x^ii:
  a_s:=a_s+a||ii*x^ii:
end do:
print(a_s);
print(w_s);
for ii in H,q,Y do
  tmp[dune_series][ii]:=collect(subs(w(x)=w_s,a(x)=a_s,tmp[e][ii]),x):
end do:
for ii from 0 to 10 do
  set_||ii:=[coeff(tmp[dune_series][Y],x,ii),coeff(tmp[dune_series][H],x,ii),coeff(tmp[dune_series][q],x,ii)]:
end do:

#a0:=0:
a0:=0:
w0:=wo:
#l:=0:
#M:=0:
#L:=0:
#R:=0:
knowns:={l,wo,R,L,M}:
for ii from 0 to 10 do
  unknowns_||ii:=indets(set_||ii) minus knowns:
end do:
print(solve(set_0,unknowns_0));
sol_0:=solve(set_0,unknowns_0)[1]:
#print(allvalues(solve(set_0,unknowns_0))[2]);
#sol_0:=allvalues(solve(set_0,unknowns_0))[2]:
print(sol_0);
#fin();
for jj from 1 to 10 do
  set_||jj:=simplify(eval(subs(op(sol_0),set_||jj))):
  unknowns_||jj:=indets(set_||jj,name) minus knowns:
end do:

#fin();

omnisol:={va0=a0,va1=a1,va2=a2,va3=a3,va4=a4,va5=a5,vw0=w0,vw1=w1,vw2=w2,vw3=w3,vw4=w4,vw5=w5}:
update_omnisol:=proc(sol::set)::NULL:
  global omnisol:
  local tmp,ii:
  tmp:={}:
  for ii in omnisol do
    tmp:=tmp union {simplify(subs(op(sol),ii))}:
  end do:
  omnisol:=tmp:
  print(omnisol);
  return NULL:
end proc:

update_omnisol(sol_0):
toten:=10:
for ii from 1 to toten do
  dit({fyellow,underline},"solving term %d",ii);
  dit({italic,fmagenta},"Here are the unknowns:");
  print(unknowns_||ii);
  dit({italic,fgreen},"Here are all the solutions:");
  allsols:=solve(set_||ii,unknowns_||ii):
  print(allsols);
  dit({italic,fred},"Here are the chosen solutions:");
  sol_||ii:=allsols:
  print(sol_||ii);
  for jj from ii+1 to toten do
    set_||jj:=simplify(subs(op(sol_||ii),set_||jj)):
    unknowns_||jj:=indets(set_||jj,name) minus knowns:
  end do:
  update_omnisol(sol_||ii):
  dit({italic,fcyan},"Here is the test:");
  for kk in set_||ii do
    print(simplify(eval(subs(op(sol_||ii),kk))));
  end do:
end do:
#*)

#fin();
dit({},"exporting results to python");

eprint(tmp,f,"what we are looking at");

tmp[f][H]:=tmp[f][H]-diff(J(x),x):

deqns:={tmp[f][H],tmp[f][q],tmp[f][Y]}:
inits:={}:
vars:={w(x),a(x),J(x)}:
ivar:=x:
first_order:=convertsys(deqns,inits,vars,ivar):
dit({italic,fred},"Variable definitions:");
print(first_order[2]);

tmp1:=map(rhs,first_order[1]):
str:=py_syntax(tmp1);
fprintf("late_time_first_order_system.py",cat("def theory(x,Y):\n","    return ",str)):

replace_term:=proc(term)
  return subs(op(omnisol),term):
end proc:

str:=py_syntax(map(replace_term,[J0,va0+va1*x+va2*x^2+va3*x^3,va1+2*va2*x+3*va3*x^2+4*va4*x^3,vw0+vw1*x+vw2*x^2+vw3*x^3,vw1+2*vw2*x+3*vw3*x^2+4*vw4*x^3]));
fprintf("late_time_boundary.py",cat("Y0=",str)):
fin();



fin();














#	again, a successful attempt to incorporate ad hoc cosmological constant

read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:

dit({},"First we would like to check (88) for Class 3C:");

for ii in H,q,X,Y do
  tmp[a][ii]:=simplify(subs(A=0,S3=0,S2=-1/3,S1=-1/3,U2=-4/3,U1=-K*L/3,Q(t)=x(t)/sqrt(K),tmp[Friedmann_form][ii])):
end do:

K:=1:

eprint(tmp,a,"eoms");

simplify(subs(isolate(tmp[a][Y],diff(U(t),t)),tmp[a][q]));

dit({},"now reduce order");

for ii in H,q,X,Y do
  tmp[b][ii]:=simplify(subs(x(t)=sqrt(y(t)),O_d(t)=K*Rm(t)/(3*H(t)^2),O_r(t)=K*Rr(t)/(3*H(t)^2),O_L(t)=l/(3*H(t)^2),tmp[a][ii])):
end do:

tmp[b][Y]:=simplify(y(t)^(3/2)*tmp[b][Y]):
tmp[b][H]:=simplify(y(t)*H(t)^2*tmp[b][H]):
tmp[b][q]:=simplify(y(t)*H(t)^2*tmp[b][q]):

for ii in H,q,Y,X do
  if ii<>X then
    tmp[c][ii]:=simplify(subs(isolate(tmp[b][X],diff(y(t),t)),tmp[b][ii])):
    tmp[c][ii]:=simplify(subs(isolate(tmp[b][X],diff(y(t),t)),tmp[c][ii])):
  else
    tmp[c][X]:=tmp[b][X]:
  end if:
end do:

eprint(tmp,c,"the neatest first-order form available");

dit({},"maple attempt at the first order system");

tmp[c][J]:=diff(J(t),t)=lhs(tmp[c][H])-rhs(tmp[c][H]):
tmp[c][Rm]:=diff(Rm(t),t)=-3*Rm(t)*H(t):
tmp[c][Rr]:=diff(Rr(t),t)=-4*Rr(t)*H(t):

deqns:={tmp[c][X],tmp[c][Y],tmp[c][q],tmp[c][Rm],tmp[c][Rr],tmp[c][J]}:
inits:={}:
vars:={U(t),y(t),H(t),Rm(t),Rr(t),J(t)}:
ivar:=t:
first_order:=convertsys(deqns,inits,vars,ivar):
dit({},"Here are the variable definitions:");
print(first_order[2]);
tmp1:=map(rhs,first_order[1]):
str:=py_syntax(tmp1):
fprintf("bounce.py",cat("def theory(Y,t):\n","    return ",str)):

dit({},"The integral constraint:");

print(tmp[c][H]);
sol:=solve(tmp[c][H],y(t));

str:=py_syntax(subs(U(t)=U,H(t)=H,Rm(t)=Rm,Rr(t)=Rr,sol[1])):
fprintf("constraint_1.py",cat("def constraint1(H,U,L,l,Rm,Rr):\n","    return ",str)):
str:=py_syntax(subs(U(t)=U,H(t)=H,Rm(t)=Rm,Rr(t)=Rr,sol[2])):
fprintf("constraint_2.py",cat("def constraint2(H,U,L,l,Rm,Rr):\n","    return ",str)):
print(str);

fin();

#	this is the successful attempt to include matter and radiation in the inflationary formalism

read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:

dit({},"First we would like to check (88) for Class 3C:");

for ii in H,q,X,Y do
  tmp[a][ii]:=simplify(subs(A=0,S3=0,S2=-1/3,S1=-1/3,U2=-4/3,U1=-K*L/3,Q(t)=x(t)/sqrt(K),tmp[Friedmann_form][ii])):
end do:

K:=1:

eprint(tmp,a,"eoms");

simplify(subs(isolate(tmp[a][Y],diff(U(t),t)),tmp[a][q]));

dit({},"now reduce order");

for ii in H,q,X,Y do
  tmp[b][ii]:=simplify(subs(x(t)=sqrt(y(t)),O_d(t)=K*Rm(t)/(3*H(t)^2),O_r(t)=K*Rr(t)/(3*H(t)^2),O_L(t)=0,tmp[a][ii])):
end do:

tmp[b][Y]:=simplify(y(t)^(3/2)*tmp[b][Y]):
tmp[b][H]:=simplify(y(t)*H(t)^2*tmp[b][H]):
tmp[b][q]:=simplify(y(t)*H(t)^2*tmp[b][q]):

for ii in H,q,Y,X do
  if ii<>X then
    tmp[c][ii]:=simplify(subs(isolate(tmp[b][X],diff(y(t),t)),tmp[b][ii])):
    tmp[c][ii]:=simplify(subs(isolate(tmp[b][X],diff(y(t),t)),tmp[c][ii])):
  else
    tmp[c][X]:=tmp[b][X]:
  end if:
end do:

eprint(tmp,c,"the neatest first-order form available");

dit({},"maple attempt at the first order system");

tmp[c][J]:=diff(J(t),t)=lhs(tmp[c][H])-rhs(tmp[c][H]):
tmp[c][Rm]:=diff(Rm(t),t)=-3*Rm(t)*H(t):
tmp[c][Rr]:=diff(Rr(t),t)=-4*Rr(t)*H(t):

deqns:={tmp[c][X],tmp[c][Y],tmp[c][q],tmp[c][Rm],tmp[c][Rr],tmp[c][J]}:
inits:={}:
vars:={U(t),y(t),H(t),Rm(t),Rr(t),J(t)}:
ivar:=t:
first_order:=convertsys(deqns,inits,vars,ivar):
dit({},"Here are the variable definitions:");
print(first_order[2]);
tmp1:=map(rhs,first_order[1]):
str:=py_syntax(tmp1):
fprintf("bounce.py",cat("def theory(Y,t):\n","    return ",str)):

dit({},"The integral constraint:");

print(tmp[c][H]);
sol:=solve(tmp[c][H],y(t));

str:=py_syntax(subs(U(t)=U,H(t)=H,Rm(t)=Rm,Rr(t)=Rr,sol[1])):
fprintf("constraint_1.py",cat("def constraint1(H,U,L,Rm,Rr):\n","    return ",str)):
str:=py_syntax(subs(U(t)=U,H(t)=H,Rm(t)=Rm,Rr(t)=Rr,sol[2])):
fprintf("constraint_2.py",cat("def constraint2(H,U,L,Rm,Rr):\n","    return ",str)):
print(str);

fin();









#	here is the very successful attempt at matter/radiation-free early inflation with U1 or Lambda

read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:

dit({},"First we would like to check (88) for Class 3C:");

for ii in H,q,X,Y do
  tmp[a][ii]:=simplify(subs(A=0,S3=0,S2=-1/3,S1=-1/3,U2=-4/3,U1=-K*L/3,Q(t)=x(t)/sqrt(K),tmp[Friedmann_form][ii])):
end do:

K:=1:

eprint(tmp,a,"eoms");

simplify(subs(isolate(tmp[a][Y],diff(U(t),t)),tmp[a][q]));
#dsolve(simplify(subs(H(t)=H,tmp[a][Y])),U(t));

dit({},"now reduce order");

for ii in H,q,X,Y do
  tmp[b][ii]:=simplify(subs(x(t)=sqrt(y(t)),O_d(t)=0,O_r(t)=0,O_L(t)=0,tmp[a][ii])):
end do:

tmp[b][Y]:=simplify(y(t)^(3/2)*tmp[b][Y]):
tmp[b][H]:=simplify(y(t)*H(t)^2*tmp[b][H]):
tmp[b][q]:=simplify(y(t)*H(t)^2*tmp[b][q]):


for ii in H,q,Y,X do
  if ii<>X then
    tmp[c][ii]:=simplify(subs(isolate(tmp[b][X],diff(y(t),t)),tmp[b][ii])):
    tmp[c][ii]:=simplify(subs(isolate(tmp[b][X],diff(y(t),t)),tmp[c][ii])):
  else
    tmp[c][X]:=tmp[b][X]:
  end if:
end do:

eprint(tmp,c,"the neatest first-order form available");

dit({},"maple attempt at the first order system");

tmp[c][J]:=diff(J(t),t)=lhs(tmp[c][H]):

deqns:={tmp[c][X],tmp[c][Y],tmp[c][q],tmp[c][J]}:
inits:={}:
vars:={U(t),y(t),H(t),J(t)}:
ivar:=t:
first_order:=convertsys(deqns,inits,vars,ivar):
dit({},"Here are the variable definitions:");
print(first_order[2]);
tmp1:=map(rhs,first_order[1]):
str:=py_syntax(tmp1):
fprintf("bounce.py",cat("def theory(Y,t):\n","    return ",str)):

dit({},"it would be nice to pin the Q value...");

print(tmp[c][X]);

dit({},"The integral constraint:");

print(tmp[c][H]);
sol:=solve(tmp[c][H],y(t));

str:=py_syntax(subs(U(t)=U,H(t)=H,sol[1])):
fprintf("constraint_1.py",cat("def constraint1(H,U,L):\n","    return ",str)):
str:=py_syntax(subs(U(t)=U,H(t)=H,sol[2])):
fprintf("constraint_2.py",cat("def constraint2(H,U,L):\n","    return ",str)):

fin();








