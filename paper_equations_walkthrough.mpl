with(CLIo):
read `tools.mpl`:
with(tools):
read `theory_tools.mpl`:

dit({},"Next we would like to check (43):");

convert_parameters_print_equations(W,M);

dit({},"Next we would like to check (50):");

convert_parameters_print_equations(NY,W);
convert_parameters_print_equations(Y,W);

dit({},"Next we would like to check (68):");

alias(a=a(t),b=b(t),X=X(t),Y=Y(t)):
read `reduced_PGT_Lagrangian.mpl`:
full_pgt_lag:=%:
alias(a='a',b='b',X='X',Y='Y'):
full_pgt_lag:=simplify(subs(A1=0,A2=0,A3=0,A4=0,A5=0,A6=0,B1=0,B2=0,B3=0,full_pgt_lag)):
alias(R=a(t),S=b(t),X=X(t),Y=Y(t)):
print(full_pgt_lag);
dit({},"coeff of A in reduced lagrangian:");
print(coeff(full_pgt_lag,A));
alias(R='R',S='S',X='X',Y='Y'):

dit({},"Next we would like to check (78) and (80):");

convert_parameters_print_equations(cW2,W);
convert_parameters_print_equations(cW1,W);

read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:

dit({},"Next (84a)-(84d) full equations:");
read `forms/conformalXY_form.mpl`:
tmp[conformalXY_form]:=%:
for ii in X,Y,a,b do
  tmp[conformalXY_form2][ii]:=simplify(expand(lhs(tmp[conformalXY_form][ii])),size):
end do:

tmp[conformalXY_form2][X]:=-tmp[conformalXY_form2][X]:
tmp[conformalXY_form2][Y]:=-tmp[conformalXY_form2][Y]:

alias(X=X(u),Y=Y(u),R=a(u)):
eprint(tmp,conformalXY_form2,"conformal equations:");
dit({},"the coefficient of S3 in the a equation:");
print(simplify(coeff(tmp[conformalXY_form2][a],S3),size)):
dit({},"the coefficient of S2 in the a equation:");
print(simplify(coeff(tmp[conformalXY_form2][a],S2),size)):
alias(X='X',Y='Y',R='R'):

dit({},"First we would like to check (88) for Class 3C:");

for ii in H,q,X,Y do
  tmp[Class_3C][ii]:=simplify(subs(A=0,U1=0,S3=0,tmp[Friedmann_form][ii])):
end do:

U_sol:=isolate(tmp[Class_3C][X],U(t)):

print(U_sol);

dit({},"Next we would like to find the modified gravitational dimensionless densities in (94)");

`latex/special_names`["S1"]:="\\sigma_1":
`latex/special_names`["S2"]:="\\sigma_2":
`latex/special_names`["S3"]:="\\sigma_3":
`latex/special_names`["U1"]:="\\upsilon_1":
`latex/special_names`["U2"]:="\\upsilon_2":
`latex/special_names`["K"]:="\\kappa":
`latex/special_names`["dQ"]:="\\partial_{t}Q":

H_eqn:=simplify(subs(U_sol,tmp[Class_3C][H])):
H_eqn:=simplify(rhs(H_eqn)-lhs(H_eqn)):
sum_new_Omegas:=H_eqn-O_d(t)-O_r(t)-O_L(t):
sum_new_Omegas:=subs(diff(Q(t),t)=H(t)*x,sum_new_Omegas):
sum_new_Omegas:=collect(sum_new_Omegas,x,simplify):
sum_new_Omegas:=subs(x=dQ/H(t),sum_new_Omegas):
sum_new_Omegas:=subs(Q(t)=Q,sum_new_Omegas):
sum_new_Omegas:=subs(H(t)=H,sum_new_Omegas):
print(sum_new_Omegas);
print(op(1,sum_new_Omegas));
print(op(2,sum_new_Omegas));
print(op(3,sum_new_Omegas));

#sum_new_Omegas=algsubs(Q(t)=Q,sum_new_Omegas):
latex(op(1,sum_new_Omegas),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/modified_gravitational_densities_1.tex"):
latex(op(2,sum_new_Omegas),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/modified_gravitational_densities_2.tex"):
latex(op(3,sum_new_Omegas),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/modified_gravitational_densities_3.tex"):


dit({},"Next we would like to obtain the f_i coefficients in (96)");

Q_eqn:=subs(U_sol,tmp[Class_3C][Y]):
Q_eqn2:=simplify(numer(normal(lhs(Q_eqn)-rhs(Q_eqn)))):
print(Q_eqn2);
inert_auxiliary:=simplify(subs(diff(diff(Q(t),t),t)=Q(t)*x^5,Q_eqn2)):
inert_auxiliary:=simplify(subs(diff(Q(t),t)=Q(t)*y^2,inert_auxiliary)):
inert_auxiliary:=simplify(subs(diff(H(t),t)=x^2,inert_auxiliary)):
inert_auxiliary:=simplify(subs(H(t)=x^3/y^2,inert_auxiliary)):
inert_auxiliary:=collect(inert_auxiliary,y,simplify):
inert_auxiliary:=algsubs(y^4=x^4,inert_auxiliary):
inert_auxiliary:=algsubs(y^(-4)=x^(-5),inert_auxiliary):
inert_auxiliary:=collect(inert_auxiliary,x,simplify):
print(inert_auxiliary);


(*
dit({},"Let's check that this gives the original form:");

Q_eqn2_check:=algsubs(x^5=diff(diff(Q(t),t),t)/Q(t),inert_auxiliary):
Q_eqn2_check:=algsubs(x^4=(diff(Q(t),t)/Q(t))^2,Q_eqn2_check):
Q_eqn2_check:=algsubs(x^3=H(t)*diff(Q(t),t)/Q(t),Q_eqn2_check):
Q_eqn2_check:=algsubs(x^2=diff(H(t),t),Q_eqn2_check):
Q_eqn2_check:=algsubs(x=H(t)^2,Q_eqn2_check):
print(simplify(Q_eqn2-Q_eqn2_check));
*)

inert_auxiliary:=algsubs(Q(t)=Q,inert_auxiliary):


(*
print(eval(`latex/special_names`));
`latex/Q`:=proc(t)
  sprintf("Q");
end proc:
`latex/U`:=proc(t)
  sprintf("U");
end proc:
`latex/H`:=proc(t)
  sprintf("H");
end proc:
*)

latex(simplify(coeff(inert_auxiliary,x^5)),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/f1.tex"):
latex(simplify(coeff(inert_auxiliary,x^4)),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/f2.tex"):
latex(simplify(coeff(inert_auxiliary,x^3)),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/f3.tex"):
latex(simplify(coeff(inert_auxiliary,x^2)),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/f4.tex"):
latex(simplify(coeff(inert_auxiliary,x^1)),"/home/williamb/Documents/physics/papers/paper-2/paper/aps/f5.tex"):



dit({},"Next we would like to check (97) for Class 3C:");

Q_eqn:=numer(normal(simplify(subs(S1=S2*S,H(t)=2/(3*(1+w)*t),Q(t)=sqrt(Q2),Q(t)*Q_eqn2)))):
Q_eqn:=collect(Q_eqn/Q2,Q2):
Q_eqn:=simplify(subs(Q2=KQ2dU2*U2/K,Q_eqn)):
posroot:=solve(Q_eqn,KQ2dU2)[1]:
print(posroot);

dit({},"Next we would like to check the infamous (103) for Class 3C*:");

posroot_matter:=simplify(subs(w=0,S=1,posroot));

H_eqn:=simplify(subs(U_sol,tmp[Class_3C][H])):
H_eqn:=simplify(subs(diff(Q(t),t)=0,H_eqn)):
H_eqn:=simplify(subs(S1=S2,H(t)=2/(3*(1+w)*t),Q(t)=sqrt(Q2),H_eqn)):
H_eqn:=simplify(subs(Q2=U2/(4*K*S2),H_eqn));

dit({},"Next we would like to check the renormalised Einstein constant in Fig 3 for Class 3C:");

Q_eqn:=subs(U_sol,tmp[Class_3C][Y]):
Q_eqn2:=simplify(numer(normal(lhs(Q_eqn)-rhs(Q_eqn)))):
Q_eqn:=numer(normal(simplify(subs(S1=S2*S,H(t)=1/(2*t),Q(t)=sqrt(Q2),Q(t)*Q_eqn2)))):
Q_eqn:=collect(Q_eqn/Q2,Q2):
Q_eqn:=simplify(subs(Q2=KQ2dU2*U2/K,Q_eqn));
posroot:=solve(Q_eqn,KQ2dU2):
print(posroot);

H_eqn:=simplify(subs(U_sol,tmp[Class_3C][H])):
H_eqn:=simplify(subs(diff(Q(t),t)=0,H_eqn)):
H_eqn:=simplify(subs(S1=S2*S,H(t)=1/(2*t),Q(t)=sqrt(Q2),H_eqn));
H_eqn:=simplify(subs(Q2=U2/(4*K*S2*S),H_eqn));

dit({},"Next we would like to check (118a) and (118b) for Class 4H:");

for ii in H,q,X,Y do
  tmp[Class_4H][ii]:=simplify(subs(A=0,U1=0,S3=0,S2=0,tmp[Friedmann_form][ii])):
end do:

H_sol:=isolate(tmp[Class_4H][X],H(t)):

print(simplify(sqrt(3)*H_sol));

print(tmp[Class_4H][Y]);

print(simplify(subs(H_sol,3*H(t)-U(t))/sqrt(3)));

sum_new_rhos:=simplify(subs(H_sol,-3*H(t)^2*lhs(tmp[Class_4H][H])/(K*U2)));
sum_new_rhos:=simplify(subs(U(t)=3*H(t)-chi/R(t),lhs(tmp[Class_4H][H])));

dit({},"Next we would like to check (118a) extended to Class 4I:");

for ii in H,q,X,Y do
  tmp[Class_4I][ii]:=simplify(subs(A=0,U1=0,S3=0,U2=0,tmp[Friedmann_form][ii])):
end do:

print(isolate(tmp[Class_4I][X],H(t)));

dit({},"Next we would like to check (128) for Class 3E:");

for ii in H,q,X,Y do
  tmp[Class_3E][ii]:=simplify(subs(A=0,U2=0,S3=0,O_L(t)=0,tmp[Friedmann_form][ii])):
end do:

Q_sol:=simplify(tmp[Class_3E][H]-tmp[Class_3E][q]):

print(Q_sol);

dit({},"Next we would like to check (129) and (130) for Class 3E:");

rels:=O_d(t)=H_0^2*O_d0/(diff(a(t),t)*a(t))^2,O_r(t)=H_0^2*O_r0/(diff(a(t),t)^2*a(t)^3),H(t)=diff(a(t),t)/a(t):

Q_sol:=simplify(subs(rels,Q_sol)):

Q_sol:=isolate(Q_sol,Q(t)^2):

Qd_sol:=isolate(diff(Q_sol,t),diff(Q(t),t)):

U_sol:=isolate(tmp[Class_3E][X],U(t)):

Q_a_eqn:=simplify(subs(rels,simplify(subs(U_sol,simplify(tmp[Class_3E][H]-2*tmp[Class_3E][q]))))):
tmpp:=simplify(algsubs(Q_sol,simplify(subs(Qd_sol,Q_a_eqn)),exact)):
Q_a_eqn:=simplify(subs(Q_sol,numer(normal(lhs(tmpp)-rhs(tmpp))))):

solve(simplify(eval(subs(a(t)=(O_r0/O_d0)*(cosh(c2*t)-1),Q_a_eqn)),trig),c2);

dit({},"Finally appendix B the Heineke theory:");

read `forms/conformalXY_form.mpl`:

tmp[conformalXY_form]:=%:
for ii in X,Y,a,b do
  tmpp:=tmp[conformalXY_form][ii]:
  tmp[conformalXY_form2][ii]:=simplify(expand(subs(U1=0,U2=0,A=0,Ls=0,Rho_d=0,S1=0,S3=S2/3,tmpp)),size):
end do:

alias(X=X(u),Y=Y(u),R=a(u)):
eprint(tmp,conformalXY_form2,"conformal equations:");
alias(X='X',Y='Y',R='R'):

for ii in H,q,X,Y do
  tmp[Class_s][ii]:=simplify(subs(Q(t)=0,U(t)=0,U1=0,U2=0,A=0,O_L(t)=0,O_d(t)=0,S1=0,S3=S2/3,tmp[Friedmann_form][ii])):
end do:
eprint(tmp,Class_s,"the Friedmann form for the silly class:");

dit({},"Next we would like to check (132,133,134)");

dit({},"D");
testcase_critical_case(16);
dit({},"H");
testcase_critical_case(14);
dit({},"J");
testcase_critical_case(11);
dit({},"O");
testcase_critical_case(10);
dit({},"E");
testcase_critical_case(1);
dit({},"E extras");
testcase_critical_case(27);
testcase_critical_case(30);
testcase_critical_case(35);
dit({},"of interest...");
testcase_critical_case(15);
testcase_critical_case(12);

dit({},"Finally appendix B all the literature constraints:");

testcase_literature(Minkevich);
testcase_literature(Zhang);
testcase_literature(SNY1);
testcase_literature(SNY2);

fin();
