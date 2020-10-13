comment({fblack,bmagenta,bold,underline},"This is hopefully just an interlude to check if we are talking about the same system as Anthony...");

#tmp:=table([]);
tmp[Anthony_raw][Hderiv]:=diff(H(t),t)=-1/324*(-216*Phi0^5*lambda*Q0(t)^2*sigma2*eta*sigma1+1296*Phi0*eta*H(t)^2*Q0(t)^4*sigma2*sigma1^2-324*Phi0*eta*H(t)^2*Q0(t)^4*sigma2^2*sigma1-36*Phi0*eta*Q0(t)^4*sigma2^3*U0(t)^2-144*Phi0*eta*Q0(t)^4*U0(t)^2*sigma1^3+6*Phi0^3*eta^2*Q0(t)^2*sigma2*U0(t)^2*sigma1+216*Phi0*eta*Q0(t)^4*sigma2^3*U0(t)*H(t)+144*Phi0*eta*Q0(t)^4*U0(t)^2*sigma1^2*sigma2+36*Phi0*eta*Q0(t)^4*U0(t)^2*sigma1*sigma2^2+864*Phi0*eta*Q0(t)^4*U0(t)*H(t)*sigma1^3+10368*rho(t)*Pi*Q0(t)^4*sigma1^3+5184*Phi0^3*lambda*Q0(t)^4*sigma1^3-31104*P(t)*Pi*Q0(t)^4*sigma1^3-324*Phi0*eta*H(t)^2*Q0(t)^4*sigma2^3-2592*rho(t)*Pi*Q0(t)^4*sigma1*sigma2^2-1296*Phi0^3*lambda*Q0(t)^4*sigma1*sigma2^2+7776*P(t)*Pi*Q0(t)^4*sigma1*sigma2^2+1296*Phi0*eta*H(t)^2*Q0(t)^4*sigma1^3-432*rho(t)*Pi*Q0(t)^2*sigma2*eta*Phi0^2*sigma1+Phi0^5*eta^3*sigma2*U0(t)^2+1296*P(t)*Pi*Q0(t)^2*sigma2*eta*Phi0^2*sigma1-864*Phi0*eta*Q0(t)^4*U0(t)*H(t)*sigma1^2*sigma2-216*Phi0*eta*Q0(t)^4*U0(t)*H(t)*sigma1*sigma2^2)/Phi0/eta/Q0(t)^4/(2*sigma1-sigma2)/sigma1/(2*sigma1+sigma2);
tmp[Anthony_raw][Qderiv]:=diff(Q0(t),t)=-1/18*(U0(t)*Phi0^2*eta+6*U0(t)*Q0(t)^2*sigma2-18*Q0(t)^2*H(t)*sigma2+18*Q0(t)^2*H(t)*sigma1)/Q0(t)/sigma1;
tmp[Anthony_raw][Uderiv]:=diff(U0(t),t)=-1/3*(72*rho(t)*Pi-216*P(t)*Pi+9*Phi0*eta*U0(t)*H(t)+36*Phi0^3*lambda-Phi0*eta*U0(t)^2)/eta/Phi0;
tmp[Anthony_raw][constraintt]:=324*Q0(t)^4*sigma2^3*U0(t)*H(t)-3/2*sigma2*Phi0^4*eta^2*U0(t)^2+1944*sigma2*Q0(t)^4*H(t)^2*sigma1^2-648*Q0(t)^2*Phi0^4*lambda*sigma1^2-5184*rhor(t)*Pi*Q0(t)^2*sigma1^2+216*Q0(t)^4*sigma2*U0(t)^2*sigma1^2-18*Q0(t)^2*sigma2^2*eta*Phi0^2*U0(t)^2+36*Q0(t)^2*eta*Phi0^2*U0(t)^2*sigma1^2-5184*Q0(t)^2*Pi*Phi0*rho(t)*sigma1^2-1296*Q0(t)^4*H(t)*sigma2*U0(t)*sigma1^2+54*sigma2^2*Phi0^2*eta*H(t)*Q0(t)^2*U0(t)-216*Q0(t)^2*eta*Phi0^2*U0(t)*H(t)*sigma1^2-486*sigma2^3*Q0(t)^4*H(t)^2-54*Q0(t)^4*sigma2^3*U0(t)^2=0;

eprint(tmp,Anthony_raw,"So, here are Anthony's equations:");

#	let's put everything in my form except Anthony's eWGT coeffs for the matter part...

rels:=U0(t)=U(t),Q0(t)=Q(t),sigma1=-S1,sigma2=-S2,Phi0=1/(sqrt(K)),eta=3*U3,rho(t)=3*H(t)^2*O_d(t)/(8*Pi),rhor(t)=3*H(t)^2*O_r(t)/(8*Pi),lambda=3*H(t)^2*O_L(t),P(t)=0;

for ii in Hderiv,Qderiv,Uderiv,constraintt do
  tmp[fullconvert][ii]:=simplify(subs(rels,tmp[Anthony_raw][ii]));
end do;

eprint(tmp,fullconvert,"now having put in our terms...");



solve(tmp[fullconvert][Qderiv],U(t));

tmp[ants][Usol]:=U(t)=%;

for ii in Hderiv,Uderiv,constraintt do
 tmp[ants][ii]:=simplify(subs(tmp[ants][Usol],tmp[fullconvert][ii]));
end do;

eprint(tmp,ants,"we have substituted for the U solution:");
for ii in constraintt,Uderiv,Hderiv do
  tmp[all][ii]:=simplify(primpart(numer(normal(lhs(tmp[ants][ii])-rhs(tmp[ants][ii])))))=0;
end do;

for ii in constraintt,Hderiv do
  tmp[source][ii]:=tmp[all][ii];
end do;

K:=1;
eprint(tmp,all,"perhaps a clearer format:");

#W();

shuffle:=proc(form,eqs::table,version::name)::equation;
  local ii,wanted_coeffs,tmp,wanted_functions,nsol,wanted_coeffs_v,wanted_functions_v,labels,coef,buildvec,sys,sol,filler,n;
  wanted_coeffs:=proc(ii::integer)::rational;
    return sign(op(ii,form))*content(op(ii,form));
  end proc;
  wanted_functions:=proc(ii::integer)::function;
    return sign(op(ii,form))*primpart(op(ii,form));
  end proc; 
  wanted_coeffs_v:=Vector(nops(form),wanted_coeffs);
  wanted_functions_v:=Vector(nops(form),wanted_functions);
  labels:=map(op,[indices(eqs[version])]);
  for ii in op(labels) do
    tmp[zeros][ii]:=simplify(lhs(eqs[version][ii])-rhs(eqs[version][ii]));
  end do;
  coef:=proc(funct::function,target);
    return coeff(target,funct);
  end proc;
  buildvec:=proc(ii::name)::vector;
    return map(coef,wanted_functions_v,tmp[zeros][ii]);
  end proc;
  sys:=map(buildvec,labels);
  sys:=convert(sys,Matrix);
  sol:=convert(LinearAlgebra[LinearSolve](sys,wanted_coeffs_v),list); 
  comment({},"solution");
  print(sol);
  nsol:=LinearAlgebra[NullSpace](sys);
  if nsol<>{} then
    sol:=convert(nsol[1],list);
  end if;
  comment({},"null space");
  print(sol);
  filler:=[];
  n:=0;
  for ii in op(labels) do
    n:=n+1;
    filler:=[op(filler),ii=sol[n]];
  end do;
  return table([op(filler)]);
end proc;
form_H:=O_d(t)+O_r(t)+O_L(t);
H_coeffs:=shuffle(form_H,tmp,source);
form_q:=O_d(t)/2+O_r(t)-O_L(t);
q_coeffs:=shuffle(form_q,tmp,source);
initial_shuffle:=proc(form,resl,eqns)::equation;
  local sm;
  sm:=0;
  for ii in eqns do
  sm:=sm+resl[ii]*tmp[all][ii];
  end do;
  if evalb(form=0) then
    sm:=simplify(numer(normal(lhs(sm))))=0;
  else
    sm:=simplify(solve(algsubs(form=x,sm),x))=form;
  end if;
  return sm;
end proc;
tmp[shuffle][H]:=initial_shuffle(form_H,H_coeffs,{constraintt,Hderiv});
tmp[shuffle][q]:=initial_shuffle(form_q,q_coeffs,{constraintt,Hderiv});
auxiliary_coeffs:=shuffle(form_H,tmp,all);
tmp[shuffle][auxiliary]:=initial_shuffle(0,auxiliary_coeffs,{constraintt,Uderiv,Hderiv});
#tmp[shuffle][Y]:=tmp[omegas][Y];
tmp[shuffle][U]:=tmp[ants][Usol];
eprint(tmp,shuffle,"This is the very basic attempt at shuffling");

clarify:=proc(eq::equation)::equation;
  local lh;
  lh:=lhs(eq);
  expand(%);
  simplify(%);
  normal(%);
  simplify(%);
  return %=rhs(eq);
end proc;

tmp[Anthony][H]:=tmp[shuffle][H];
tmp[Anthony][q]:=tmp[shuffle][q];
tmp[Anthony][U]:=tmp[shuffle][U];
tmp[Anthony][auxiliary]:=tmp[shuffle][auxiliary];

tools[eprint2]:=proc(tab::name,eqs1::name,eqs2::name,commentry::string:="(no commentary provided to eprint!)")::NULL;
  local ii;
  comment({underline,fblack,bblue,bold},commentry);
  for ii in op(map(op,[indices(tab[eqs1])])) do
    comment({},"%s-equation in %s:\n",ii,eqs1);
    print(clarify(tab[eqs1][ii]));
    comment({},"%s-equation in %s:\n",ii,eqs2);
    print(clarify(tab[eqs2][ii]));
  end do;
  return NULL;
end proc;

K:='K';

tmp[Will][H]:=tmp[U_free][H];
tmp[Will][q]:=tmp[U_free][q];
tmp[Will][U]:=U_sol;
tmp[Will][auxiliary]:=numer(normal(lhs(tmp[U_free][Y])))=0;

tools[eprint2](tmp,Will,Anthony,"COMPARE BUDDY");
