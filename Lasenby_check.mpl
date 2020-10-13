with(CLIo):
read `tools.mpl`:
with(tools):
read `theory_tools.mpl`:

dit({fblack,bblue,bold,underline},"Goenner and Muller-Hoisen:");

read `forms/Goenner_form.mpl`:
read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:
for ii in H,q,X,Y do
  tmp[Friedmann_converted][ii]:=simplify(subs(K=1/M_P^2,tmp[Friedmann_form][ii])):
end do:
for ii in density,pressure,tor1,tor2 do
  tmp[Goenner_converted][ii]:=simplify(subs(K=1/M_P^2,tmp[Goenner_form][ii])):
end do:

rels:=h(t)=-U(t)/3,f(t)=Q(t),w=0,qtor=0,stor=0,k=-a(t)^2*H(t)^2*O_k(t),rho(t)=3*H(t)^2*O_d(t)*M_P^2,c_0=0:

for ii in tor1,tor2,density,pressure do
  tmp1:=simplify(subs(rels,tmp[Goenner_converted][ii])):
  tmp1:=simplify(subs(diff(a(t),t)=a(t)*H(t),tmp1)):
  tmp1:=simplify(subs(diff(a(t),t)=a(t)*H(t),tmp1)):
  tmp1:=simplify(subs(diff(a(t),t)=a(t)*H(t),tmp1)):
  tmp[G][ii]:=simplify(numer(normal(lhs(tmp1)-rhs(tmp1)))):
end do:


print(convert_parameters_return_equations(cW1,G));
parameter_rels:=S1=3*c_5/2+c_6/4+c_7/4-c_9/4,S2=3*c_5/2+c_6/2+c_7/2-3*c_8/4+c_9/4,S3=3*c_5/2+c_6/2+c_7/2-c_8/4-c_9/4,U1=c_1+3*c_2,U2=-c_1-3*c_3,A=2*c_4:
print(parameter_rels);
rels:=O_r(t)=0,O_L(t)=0:

for ii in X,Y,H,q do
  tmp1:=simplify(subs(parameter_rels,rels,tmp[Friedmann_converted][ii])):
  tmp[F][ii]:=simplify(numer(normal(lhs(tmp1)-rhs(tmp1)))):
end do:

M_P:=1:

F_H_Qd2:=coeff(tmp[F][H],diff(Q(t),t)^2):
F_X_Udd:=coeff(tmp[F][X],diff(diff(U(t),t),t)):
M_d_Qd2:=coeff(tmp[G][density],diff(Q(t),t)^2):
M_d_Udd:=coeff(tmp[G][density],diff(diff(U(t),t),t)):
tmp[Fs][d]:=simplify((M_d_Qd2/F_H_Qd2)*tmp[F][H]+(M_d_Udd/F_X_Udd)*tmp[F][X]):
dit({},"compare density and d equations");
#print(tmp[Fs][d]):
#print(simplify(tmp[G][density])):
print(simplify(tmp[Fs][d]-tmp[G][density]));
dit({},"compare X and tor1 equations");
F_X_Hdd:=coeff(tmp[F][X],diff(diff(H(t),t),t)):
M_tor1_Hdd:=coeff(tmp[G][tor1],diff(diff(H(t),t),t)):
tmp[Fs][X]:=simplify((M_tor1_Hdd/F_X_Hdd)*tmp[F][X]):
#print(simplify(tmp[Fs][X])):
#print(simplify(tmp[G][tor1])):
print(simplify(tmp[Fs][X]-tmp[G][tor1]));
dit({},"compare Y and tor2 equations");
F_Y_Qdd:=coeff(tmp[F][Y],diff(diff(Q(t),t),t)):
M_tor2_Qdd:=coeff(tmp[G][tor2],diff(diff(Q(t),t),t)):
tmp[Fs][Y]:=(M_tor2_Qdd/F_Y_Qdd)*tmp[F][Y]:
#print(simplify(tmp[Fs][Y])):
#print(simplify(tmp[G][tor2])):
print(simplify(tmp[Fs][Y]-tmp[G][tor2]));

dit({fblack,bblue,bold,underline},"Minkevich:");

read `forms/Minkevich_form.mpl`:
read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:
for ii in H,q,X,Y do
  tmp[Friedmann_converted][ii]:=simplify(subs(K=1/M_P^2,tmp[Friedmann_form][ii])):
end do:
for ii in density,pressure,tor1,tor2 do
  tmp[Minkevich_converted][ii]:=simplify(subs(K=1/M_P^2,tmp[Minkevich_form][ii])):
end do:

rels:=S_1(t)=U(t)/6,S_2(t)=Q(t)/2,w=0,rho(t)=3*H(t)^2*O_d(t)*M_P^2:

for ii in tor1,tor2,density,pressure do
  tmp1:=simplify(subs(parameter_rels,rels,tmp[Minkevich_converted][ii])):
  tmp[M][ii]:=numer(normal(lhs(tmp1)-rhs(tmp1))):
end do:

print(translate_cosmological_parameters(cW1,cV));
parameter_rels:=S1=q_1/4,S2=f/2+q_2/2,S3=f/2,U1=(1/2)*b,U2=(1/4)*a,A=2*f_0:
print(parameter_rels);
rels:=O_k(t)=0,O_r(t)=0,O_L(t)=0:

for ii in X,Y,H,q do
  tmp1:=simplify(subs(parameter_rels,rels,tmp[Friedmann_converted][ii])):
  tmp[F][ii]:=numer(normal(lhs(tmp1)-rhs(tmp1))):
end do:

M_P:=1:

F_H_Qd2:=coeff(tmp[F][H],diff(Q(t),t)^2):
F_X_Udd:=coeff(tmp[F][X],diff(diff(U(t),t),t)):
M_d_Qd2:=coeff(tmp[M][density],diff(Q(t),t)^2):
M_d_Udd:=coeff(tmp[M][density],diff(diff(U(t),t),t)):
tmp[Fs][d]:=simplify((M_d_Qd2/F_H_Qd2)*tmp[F][H]+(M_d_Udd/F_X_Udd)*tmp[F][X]):
dit({},"compare density and d equations");
#print(tmp[Fs][d]):
#print(simplify(tmp[M][density])):
print(simplify(tmp[Fs][d]-tmp[M][density]));
dit({},"compare X and tor1 equations");
F_X_Hdd:=coeff(tmp[F][X],diff(diff(H(t),t),t)):
M_tor1_Hdd:=coeff(tmp[M][tor1],diff(diff(H(t),t),t)):
tmp[Fs][X]:=(M_tor1_Hdd/F_X_Hdd)*tmp[F][X]:
#print(simplify(tmp[Fs][X])):
#print(simplify(tmp[M][tor1])):
print(simplify(tmp[Fs][X]-tmp[M][tor1]));
dit({},"compare Y and tor2 equations");
F_Y_Qdd:=coeff(tmp[F][Y],diff(diff(Q(t),t),t)):
M_tor2_Qdd:=coeff(tmp[M][tor2],diff(diff(Q(t),t),t)):
tmp[Fs][Y]:=(M_tor2_Qdd/F_Y_Qdd)*tmp[F][Y]:
#print(simplify(tmp[Fs][Y])):
#print(simplify(tmp[M][tor2])):
print(simplify(tmp[Fs][Y]-tmp[M][tor2]));

dit({fblack,bblue,bold,underline},"Zhang:");
read `forms/Zhang_form.mpl`:
read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:
for ii in H,q,X,Y do
  tmp[Friedmann_converted][ii]:=simplify(subs(K=1/M_P^2,tmp[Friedmann_form][ii])):
end do:
for ii in density,pressure,tor1,tor2 do
  tmp[Zhang_converted][ii]:=simplify(subs(K=1/M_P^2,tmp[Zhang_form][ii])):
end do:

rels:=h(t)=U(t)/3,f(t)=Q(t)/2,w=0,rho(t)=3*H(t)^2*O_d(t)*M_P^2:

for ii in tor1,tor2,density,pressure do
  tmp1:=simplify(subs(rels,tmp[Zhang_converted][ii])):
  tmp[Z][ii]:=numer(normal(lhs(tmp1)-rhs(tmp1))):
end do:

print(translate_cosmological_parameters(cW1,cZ));
parameter_rels:=S1=-B_1/4-B_2/8+B_0/4,S2=B_0/4+B_1/2,S3=B_0/4,U1=A_0,U2=-(1/2)*A_1,A=alpha:
print(parameter_rels);
rels:=O_k(t)=0,O_r(t)=0,O_L(t)=0:

for ii in X,Y,H,q do
  tmp1:=simplify(subs(parameter_rels,rels,tmp[Friedmann_converted][ii])):
  tmp[F][ii]:=numer(normal(lhs(tmp1)-rhs(tmp1))):
end do:

M_P:=1:

F_H_Qd2:=coeff(tmp[F][H],diff(Q(t),t)^2):
F_X_Udd:=coeff(tmp[F][X],diff(diff(U(t),t),t)):
M_d_Qd2:=coeff(tmp[Z][density],diff(Q(t),t)^2):
M_d_Udd:=coeff(tmp[Z][density],diff(diff(U(t),t),t)):
tmp[Fs][d]:=simplify((M_d_Qd2/F_H_Qd2)*tmp[F][H]+(M_d_Udd/F_X_Udd)*tmp[F][X]):
dit({},"compare density and d equations");
#print(tmp[Fs][d]):
#print(simplify(tmp[Z][density])):
print(simplify(tmp[Fs][d]-tmp[Z][density]));
dit({},"compare X and tor1 equations");
F_X_Hdd:=coeff(tmp[F][X],diff(diff(H(t),t),t)):
M_tor1_Hdd:=coeff(tmp[Z][tor1],diff(diff(H(t),t),t)):
tmp[Fs][X]:=(M_tor1_Hdd/F_X_Hdd)*tmp[F][X]:
#print(simplify(tmp[Fs][X])):
#print(simplify(tmp[Z][tor1])):
print(simplify(tmp[Fs][X]-tmp[Z][tor1]));
dit({},"compare Y and tor2 equations");
F_Y_Qdd:=coeff(tmp[F][Y],diff(diff(Q(t),t),t)):
M_tor2_Qdd:=coeff(tmp[Z][tor2],diff(diff(Q(t),t),t)):
tmp[Fs][Y]:=(M_tor2_Qdd/F_Y_Qdd)*tmp[F][Y]:
#print(simplify(tmp[Fs][Y])):
#print(simplify(tmp[Z][tor2])):
print(simplify(tmp[Fs][Y]-tmp[Z][tor2]));

dit({fblack,bblue,bold,underline},"Lasenby");

M_P:='M_P':

read `forms/Lasenby_form_3.mpl`:
read `forms/Friedmann_form.mpl`:
tmp[Friedmann_form]:=%:
for ii in H,q,X,Y do
  tmp[Friedmann_converted][ii]:=simplify(subs(K=1/M_P^2,A=0,tmp[Friedmann_form][ii])):
end do:

for ii in O1,O2,h1,h2 do
  tmp1:=simplify(algsubs(exp(-2*Int(H(t),t))=R_0^2/R(t)^2,tmp[Lasenby_form][ii])):
  tmp[Lasenby_no_exponential][ii]:=simplify(algsubs(exp(-4*Int(H(t),t))=R_0^4/R(t)^4,tmp1)):
end do:

rels:=U0(t)=U(t),Q0(t)=Q(t),sigma1=-4*S1,sigma2=-4*S2,sigma3=-4*S3,Phi0=M_P,eta=6*U2,chi=-U1,rhor(t)=3*M_P^2*H(t)^2*O_r(t)/(8*Pi),rho(t)=3*M_P*H(t)^2*O_d(t)/(8*Pi),lambda=3*H(t)^2*O_L(t)/(M_P^2),P(t)=0,k=-R(t)^2*H(t)^2*O_k(t)/R_0^2:
print(rels);

for ii in O1,O2,h1,h2 do
  tmp1:=simplify(subs(rels,tmp[Lasenby_no_exponential][ii])):
  tmp[Lasenby_converted][ii]:=numer(normal(lhs(tmp1)-rhs(tmp1))):
end do:

dit({},"comparing X and O1 equations:");
tmp1:=simplify((4/(-216))*(tmp[Lasenby_converted][O1]=0)):
tmp2:=simplify(tmp[Friedmann_converted][X]):
print(simplify(tmp1-tmp2));

dit({},"comparing Y and O2 equations:");
tmp1:=simplify((-1/18)*(tmp[Lasenby_converted][O2]=0)):
tmp2:=simplify(tmp[Friedmann_converted][Y]):
print(simplify(tmp1-tmp2));

dit({},"coefficients of O_r,O_d,O_L:");
for ii in indices(tmp[Lasenby_converted]) do
  #dit({},"%s-equation:",convert(ii[1],string)):
  cr:=simplify(coeff(tmp[Lasenby_converted][ii[1]],O_r(t))):
  cd:=simplify(coeff(tmp[Lasenby_converted][ii[1]],O_d(t))):
  cL:=simplify(coeff(tmp[Lasenby_converted][ii[1]],O_L(t))):
  #print(cr,cd,cL):
  cr:=0:
  cd:=0:
  cL:=0:
end do:

dit({},"comparing H and h1 equations:");
tmp1:=simplify((1)*(tmp[Lasenby_converted][h1])):
tmp2:=simplify(numer(normal(lhs(tmp[Friedmann_converted][H])-rhs(tmp[Friedmann_converted][H])))):
print(simplify(tmp1-tmp2));

dit({},"comparing q and h1-h2 equations:");
tmp1:=simplify((1/2)*(tmp[Lasenby_converted][h1]-tmp[Lasenby_converted][h2])):
tmp2:=simplify(numer(normal(lhs(tmp[Friedmann_converted][q])-rhs(tmp[Friedmann_converted][q])))):
print(simplify(tmp1-tmp2));

fin();

