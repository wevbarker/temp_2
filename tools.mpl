#	a collection of tools for equations of motion in cosmology
read `tensor.mpl`:
with(StringTools);

tools[comment]:=proc(opts::set,message::string,n1:=NULL,n2:=NULL,n3:=NULL,n4:=NULL,n5:=NULL)::NULL;
  local tmp1,tmp2,transl,possible,violators,accepted,ii;
  transl:=table([none=0,bold=1,faint=2,italic=3,underline=4,sblink=5,fblink=6,conceal=8,cross=9,dunderline=21,fblack=30,fred=31,fgreen=32,fyellow=33,fblue=34,fmagenta=35,fcyan=36,fwhite=37,bblack=40,bred=41,bgreen=42,byellow=43,bblue=44,bmagenta=45,bcyan=46,bwhite=47,frame=51,encircle=52,overline=53]);
  possible:=map(op,{indices(transl)});
  if evalb(opts={}) then
    tmp2:="\e[1;34m";
  else
    if evalb(opts subset possible) then
      tmp2:="";
      for ii from 1 to numelems(opts) do
	tmp2:=cat(tmp2,convert(transl[opts[ii]],string),";");
      end do;
      tmp2:=cat("\e[",substring(tmp2,1..-2),"m");
    else
      violators:=convert(opts minus possible,string);
      accepted:=convert(possible,string);
      error cat("some of your formatting options are wrong, namely, ",violators," which are not in ",accepted);
    end if;
  end if;
  tmp1:=cat(tmp2,message,"\e[00m\n\n");
  printf(tmp1,n1,n2,n3,n4,n5);
  return NULL;
end proc;

tools[W]:=proc()::NULL;
  comment({byellow,fblack,bold,underline},"CALCULATION FINISHED");
  `quit`(12);
end proc;

convert_recursive:=proc(eq::equation,relations::set,old::function:=NULL)::equation;
  description "takes an equation and and a set of relations and repetatively substitutes until there is no change";
  local unfinished,current,previous;
  current:=eq;
  unfinished:=true;
  while unfinished do
   previous:=current;
   current:=expand(subs(op(relations),current));
   unfinished:=evalb(simplify(previous)<>simplify(current));
  end do;
  simplify(current);
  if evalb(old=NULL) then
    return simplify(current);
  else
    return simplify(content(lhs(%),old)=0);
  end if;
end proc;

tools[simplify_zero]:=proc(eq)
  description "takes an expression or equation and returns the relevant numerator which must be set to zero";
  local expression;
  if type(eq,equation) then
    expression:=lhs(eq)-rhs(eq);
  else
    expression:=eq;
  end if;
  return simplify(numer(normal(expression)));
end proc;

tools[convert_variables]:=proc(eq::equation,relations::set,old::function:=NULL)::equation;
  description "repeatedly applies a substitution until nothing changes";
  local unfinished,current,previous;
  current:=eq;
  unfinished:=true;
  while unfinished do
   previous:=current;
   current:=expand(subs(op(relations),current));
   unfinished:=evalb(simplify(previous)<>simplify(current));
  end do;
  simplify(current);
  if evalb(old=NULL) then
    return simplify(current);
  else
    return simplify(content(lhs(%),old)=0);
  end if;
end proc;

tools[print_for_paper]:=proc()::NULL;
  description "prints the gauge fields for the paper";
  local spherical_index,tmp,ii,jj,kk;
  comment({},"the b^a_mu fields");
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      comment({},"b^%d_%d:",ii,jj);
      print(grade(gammau||ii&.gd||jj,0)[1]);
    end do;
  end do;
  comment({},"the A^ij_mu");
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      for kk from 0 to 3 do
	if kk>jj then
	  comment({},"A^%d%d_%d",jj,kk,ii);
	  comment({fgreen},"geometric algebra:");
	  print(ds(grade((gammau||jj&^gammau||kk)&.Om||ii,0))[1]);
	  comment({fgreen},"tensor calculus:");
	  spherical_index:=index_translate(ii);
	  print(Au||jj||u||kk||d||spherical_index);
	end if;
      end do;
    end do;
  end do;
  comment({},"the calT^i_jk");
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      for kk from 0 to 3 do
	if kk>jj then
	  comment({},"calT^%d_%d%d",ii,kk,kk);
	  comment({fgreen},"geometric algebra:");
	  print(ds(grade((gammad||jj&^gammad||kk)&.calT_gammau_||ii,0))[1]);
	  comment({fgreen},"tensor calculus:");
	  spherical_index:=index_translate(ii);
	  print(calTu||ii||d||jj||d||kk);
	end if;
      end do;
    end do;
  end do;	  
  comment({},"the two quadratic invariants from the torsion perspective");
  comment({},"beta_1");
  print(beta_1);
  comment({},"beta_2");
  print(beta_2);
  comment({},"beta_3");
  print(beta_3);
  comment({},"now the general combination in Mike's coords");
  print(simplify(B1*beta_1+B2*beta_2+B3*beta_3));
  comment({},"now the general combination in Mike's coords v2");
  tmp:=simplify((B2/2-B1/2)*beta_1+(-B2)*beta_2+B3*beta_3);
  tmp:=collect(tmp,{B1,B2,B3},simplify);
  print(%);
end proc;
 
tools[Euler_Lagrange]:=proc(lagrangian,field)::expression;
  local tmp;
  tmp:=simplify(diff(diff(lagrangian,diff(field,t)),t)-diff(lagrangian,field));
  return tmp;
end proc;

tools[generalised_Euler_Lagrange]:=proc(lagrangian,field)::expression;
  local tmp;
  tmp:=simplify(-diff(diff(diff(lagrangian,diff(diff(field,t),t)),t),t)+diff(diff(lagrangian,diff(field,t)),t)-diff(lagrangian,field));
  return tmp;
end proc;

tools[prepare_reduced_eWGT_Lagrangian]:=proc()::expression;
  local ii,L_R,L_R2,L_T2,L_N,L_L,L_r,L_d;
  comment({bblue,bold,underline},"Forming contributions to the eWGT Lagrangian");
  #	Einstein-Hilbert term
  comment({},"Ricci scalar:");
  L_R:=-(1/2)*A*einstein_hilbert_;
  simplify(b^4*a^3*Phi^2*%);
  print(%);
  #     quadratic Riemann
  L_R2:=0;
  for ii from 1 to 6 do
    comment({},"Quadratic Riemann cA%d",ii);
    L_R2:=L_R2+A||ii*alpha_||ii||_;
    simplify(b^4*a^3*A||ii*alpha_||ii||_);
    print(%);
  od;
  #     quadratic torsion
  L_T2:=0;
  for ii from 1 to 3 do
    comment({},"Quadratic torsion cB%d",ii);
    L_T2:=L_T2+B||ii*beta_dagger_||ii||_;
    simplify(b^4*a^3*Phi^2*B||ii*beta_dagger_||ii||_);
    print(%);
  od;
  #     compensator kinetic
  comment({},"Compensator kinetic term:");
  L_N:=-N*nu_/2;
  simplify(-b^4*a^3*N*nu_/2);
  print(%);
  #     compensator interaction
  comment({},"Compensator interaction term:");
  L_L:=-Ls*K;
  simplify(-b^4*a^3*Phi^4*Ls*K);
  #     radiation 
  L_r:=-Rho_r/(b*a)^4;
  #	dust
  L_d:=-Rho_d/(b*a)^3;
  return b^4*a^3*(Phi^2*L_R+L_R2+Phi^2*L_T2+Phi^4*L_L+L_N+L_r+Phi*L_d);
end proc;

tools[prepare_reduced_PGT_Lagrangian]:=proc()::expression;
  local ii,L_R,L_R2,L_T2,L_ds,L_r,L_Ls;
  comment({bblue,bold,underline},"Forming contributions to the PGT Lagrangian");
  #	Einstein-Hilbert term
  comment({},"Ricci scalar:");
  L_R:=-(1/2)*A*einstein_hilbert_;
  simplify(b^4*a^3*%/K);
  print(%);
  #     quadratic Riemann
  L_R2:=0;
  for ii from 1 to 6 do
    comment({},"Quadratic Riemann cA%d",ii);
    L_R2:=L_R2+A||ii*alpha_||ii||_;
    simplify(b^4*a^3*A||ii*alpha_||ii||_);
    print(%);
  od;
  #     quadratic torsion
  L_T2:=0;
  for ii from 1 to 3 do
    comment({},"Quadratic torsion cB%d",ii);
    L_T2:=L_T2+B||ii*beta_||ii||_;
    simplify(b^4*a^3*B||ii*beta_||ii||_/K);
    print(%);
  od;
  #	cosmological constant
  comment({},"Cosmological constant:");
  L_Ls:=-Ls;
  simplify(-b^4*a^3*Ls/K);
  print(%);
  #     radiation 
  L_r:=-Rho_r/(a*b)^4;
  #	dust
  L_ds:=-Rho_d/((a*b)^3*sqrt(K));
  return b^4*a^3*(L_R/K+L_R2+L_T2/K+L_r+L_ds+L_Ls/K);
end proc;

tools[prepare_equations_of_motion]:=proc()::NULL;
  global tmp;
  comment({},"preparing equations of motion...");
  tmp:=[tools[Euler_Lagrange](reduced_eWGT_Lagrangian,X)=0,tools[Euler_Lagrange](reduced_eWGT_Lagrangian,Y)=0,tools[Euler_Lagrange](reduced_eWGT_Lagrangian,a)=0,tools[Euler_Lagrange](reduced_eWGT_Lagrangian,Phi)=0];
  writeto("eWGT_equations_of_motion.mpl");
  lprint(tmp); 
  print(`;`);
  writeto(terminal);
  tmp:=[tools[Euler_Lagrange](reduced_PGT_Lagrangian,X)=0,tools[Euler_Lagrange](reduced_PGT_Lagrangian,Y)=0,tools[Euler_Lagrange](reduced_PGT_Lagrangian,a)=0,tools[Euler_Lagrange](reduced_PGT_Lagrangian,b)=0];
  writeto("PGT_equations_of_motion.mpl");
  lprint(tmp); 
  print(`;`);
  writeto(terminal);
  comment({},"...done");
  return NULL;
end proc;

tools[prepare_reduced_Lagrangia]:=proc()::NULL;
  global reduced_eWGT_Lagrangian,reduced_PGT_Lagrangian;
  comment({},"preparing reduced Lagrangia...");
  reduced_eWGT_Lagrangian:=tools[prepare_reduced_eWGT_Lagrangian]();
  reduced_PGT_Lagrangian:=tools[prepare_reduced_PGT_Lagrangian]();
  writeto("reduced_eWGT_Lagrangian.mpl");
  lprint(reduced_eWGT_Lagrangian);
  print(`;`);
  writeto(terminal);
  writeto("reduced_PGT_Lagrangian.mpl");
  lprint(reduced_PGT_Lagrangian);
  print(`;`);
  writeto(terminal);
  comment({},"...done");
  return NULL;
end proc;

theory:=table([gauge_theory=NULL,Ricci=NULL,quadratic_Riemann=NULL,quadratic_torsion=NULL,compensator_kinetic=NULL,compensator_interaction=NULL,radiation=NULL,dust=NULL,curvature=NULL]);

tools[refine_equations_of_motion]:=proc()::NULL;
  global equations_of_motion,lagrangian;
  comment({},"refining equations of motion...");
  equations_of_motion:=table([]);
  #	first we deal with the equations of motion...
  alias(a=a(t),b=b(t),X=X(t),Y=Y(t),Phi=Phi(t)):
  read `PGT_equations_of_motion.mpl`;
  equations_of_motion[raw]:=%;
  alias(a=a,b=b,X=X,Y=Y,Phi=Phi):
  tools[refine_PGT_equations_of_motion]();
  alias(a=a(t),b=b(t),X=X(t),Y=Y(t),Phi=Phi(t)):
  read `eWGT_equations_of_motion.mpl`;
  equations_of_motion[raw]:=%;
  alias(a=a,b=b,X=X,Y=Y,Phi=Phi):
  tools[refine_eWGT_equations_of_motion]();
  tools[compare_gauge_theories]();
  tools[print_equations_of_motion]();
  #	now we want to refine the lagrangia...
(*
  alias(a=a(t),b=b(t),X=X(t),Y=Y(t),Phi=Phi(t)):
    read `reduced_PGT_Lagrangian.mpl`;
    lagrangian:=%;
    alias(a=a,b=b,X=X,Y=Y,Phi=Phi):
    tools[refine_PGT_Lagrangian]();
  elif theory[gauge_theory]=eWGT then
  end if;
  comment({},"...done");
*)
  return NULL;
end proc;

tools[print_equations_of_motion]:=proc()::NULL;
  dit({},"printing generalised Friedmann equations to file...");
  writeto("forms/Friedmann_form.mpl");
  lprint(eval(eoms[Friedmann])); 
  print(`;`);
  writeto(terminal);
  dit({},"...printing complete!");
  dit({},"printing generalised conformal XY equations to file...");
  writeto("forms/conformalXY_form.mpl");
  lprint(eval(eoms[PGT_conformal_form])); 
  print(`;`);
  writeto(terminal);
  dit({},"...printing complete!");
end proc;

tools[eprint]:=proc(tab::name,eqs::name,commentry::string:="(no commentary provided to eprint!)")::NULL;
  local ii;
  comment({underline,fblack,bblue,bold},commentry);
  for ii in op(map(op,[indices(tab[eqs])])) do
    comment({},"%s-equation:\n",ii);
    print(tab[eqs][ii]);
  end do;
  return NULL;
end proc;

tools[laprint]:=proc(tab,commentry::string:="(no commentary provided to eprint!)")::NULL;
  comment({underline,fblack,bblue,bold},commentry);
  print(tab);
  return NULL;
end proc;

tools[compare_gauge_theories]:=proc()::NULL;
  description "compares the eoms in conformal time for PGT and eWGT";
  global eoms;
  local PGT_X,PGT_Y,PGT_b,PGT_a,eWGT_X,eWGT_Y,eWGT_Phi,eWGT_a,diff_X,diff_Y,diff_a,diff_b;
  comment({byellow,fblue},"X equation");
  comment({bold,underline,fblack,bblue},"X");
  PGT_X:=6*eoms[PGT_conformal_form][X];
  eWGT_X:=eoms[eWGT_conformal_form][X];
  comment({},"PGT");
  print(PGT_X);
  comment({},"eWGT");
  print(eWGT_X);
  comment({},"difference");
  diff_X:=simplify(lhs(PGT_X)-lhs(eWGT_X));
  print(diff_X);
  comment({bold,underline,fblack,bblue},"Y");
  PGT_Y:=eoms[PGT_conformal_form][Y];
  eWGT_Y:=eoms[eWGT_conformal_form][Y];
  comment({},"PGT");
  print(PGT_Y);
  comment({},"eWGT");
  print(eWGT_Y);
  comment({},"difference");
  diff_Y:=simplify(lhs(PGT_Y)-lhs(eWGT_Y));
  print(diff_Y);
  comment({bold,underline,fblack,bblue},"b or Phi");
  PGT_b:=(1/K)*eoms[PGT_conformal_form][b];
  eWGT_Phi:=(1/sqrt(K))*eoms[eWGT_conformal_form][Phi];
  comment({},"PGT");
  print(PGT_b);
  comment({},"eWGT");
  print(eWGT_Phi);
  comment({},"difference");
  diff_b:=simplify(lhs(PGT_b)-lhs(eWGT_Phi),radical,symbolic); 
  print(diff_b);
  comment({bold,underline,fblack,bblue},"a");
  PGT_a:=eoms[PGT_conformal_form][a];
  eWGT_a:=eoms[eWGT_conformal_form][a];
  comment({},"PGT");
  print(PGT_a);
  comment({},"eWGT");
  print(eWGT_a);
  comment({},"difference");
  diff_a:=simplify(lhs(PGT_a)-lhs(eWGT_a));
  print(diff_a);
  if diff_X=0 and diff_Y=0 and diff_b=0 and diff_a=0 then
    comment({fred,bwhite,fblink,bold,underline},"COSMIC PGT<-->eWGT CORRESPONDENCE PROVEN");
  end if;
end proc;
  
tools[refine_PGT_Lagrangian]:=proc()::NULL;
  local simplifiers2,ii,tmp,composite_gauge_transform,specific_transform,convert_variables,relations,shuffle,form_H,H_coeffs,form_q,q_coeffs,initial_shuffle,simplifiers;
  global reduced_PGT_Lagrangian;
  lagr:=reduced_PGT_Lagrangian:
  laprint(lagr,"The general reduced PGT Lagrangian density");
  algsubs(6*A1+2*A2+2*A3+A4+A6=4*S3,lagr,exact);
  algsubs(6*A1+A2+A3+A5-A6=4*S1,%,exact);
  algsubs(6*A1+2*A2+2*A3+3*A4-4*A5+A6=4*S2,%,exact);
  algsubs(2*A5-A4=2*(S3-S2),%,exact);
  algsubs(3*B3-B1=U2,%,exact);
  algsubs(B1+3*B2=U1,%,exact);
  lagr=normal(%);
  laprint(lagr,"The specific Lagrangian density");
  dit({},"printing reduced PGT Lagrangian to file...");
  writeto("forms/PGT_Lagrangian.mpl");
  lprint(eval(lagr)); 
  print(`;`);
  writeto(terminal);
  dit({},"...printing complete!");
end proc;

tools[refine_PGT_equations_of_motion]:=proc()::NULL;
  description "Processing algorithm for the PGT equations";
  local tss,ii,tmp,composite_gauge_transform,specific_transform,convert_variables,relations,shuffle,form_H,H_coeffs,form_q,q_coeffs,initial_shuffle,simplifiers;
  global eoms;
  comment({bred,fyellow},"Now refining PGT equations");
  eoms:=table([]);
  tmp:=table([]);
  #	so first we get the raw versions
  eoms[raw][X]:=equations_of_motion[raw][1];
  eoms[raw][Y]:=equations_of_motion[raw][2];
  eoms[raw][a]:=equations_of_motion[raw][3];
  eoms[raw][b]:=equations_of_motion[raw][4];
  eprint(eoms,raw,"The raw equations of motion from the Lagrangian in terms of X, Y, a, b, Phi");
  #	next we want to get rid of any numerical factors
  for ii in X,Y,a,b do
    tmp[remove_prefactors][ii]:=simplify(eoms[raw][ii]/content(numer(normal(lhs(eoms[raw][ii])))));
  end do;
  #	we want to fix the position gauge to the Friedmann gauge early on, by setting the conformal isometry parameter to 1
  for ii in X,Y,a,b do
    tmp[friedmann][ii]:=simplify(algsubs(b(t)=1,tmp[remove_prefactors][ii],exact));
  end do;
  #	we would also like to convert things into sigma notation 
  for ii in X,Y,a,b do
    print("replacing...");
    algsubs(6*A1+2*A2+2*A3+A4+A6=4*S3,tmp[friedmann][ii],exact);
    algsubs(6*A1+A2+A3+A5-A6=4*S1,%,exact);
    algsubs(6*A1+2*A2+2*A3+3*A4-4*A5+A6=4*S2,%,exact);
    algsubs(2*A5-A4=2*(S3-S2),%,exact);
    algsubs(3*B3-B1=U2,%,exact);
    algsubs(B1+3*B2=U1,%,exact);
    tmp[cosmological_parameters][ii]:=%;
  end do;
  eprint(tmp,cosmological_parameters,"The equations of motion in the friedmann position gauge with time, t, having removed any numerical prefactors and converted to Will's cosmic coefficients");
  #	now for the paper, we want to show that the cosmic theory parameters are well motivated
  for ii in X,Y,a,b do
    tss:=tmp[cosmological_parameters][ii];
    tss:=simplify(subs(a(t)=a(t)*1,X(t)=Xw(t),Y(t)=Yw(t),%));
    tss:=dchange(t=f(u),tss,{u},known=f);
    tss:=evala(simplify(subs(diff(f(u),u)=a(u),tss)));
    tss:=simplify(subs(Xw(u)=X(u),Yw(u)=Y(u),tss));
    eoms[PGT_conformal_form][ii]:=simplify(numer(normal(lhs(tss))))=0;
  end do;
  eprint(eoms,PGT_conformal_form,"The conformal version of everything");
  simplifiers:={A,S1,S2,S3,U1,U2,N,L,Rho_r,Rho_d,k} minus {theory[Ricci]} minus {theory[quadratic_Riemann]} minus {theory[quadratic_torsion]} minus {theory[compensator_kinetic]} minus {theory[compensator_interaction]} minus {theory[radiation]} minus {theory[dust]} minus {theory[curvature]};
  for ii in X,Y,a,b do
    simplify(tmp[cosmological_parameters][ii],simplifiers);
    tmp[theory_constrained][ii]:=%;
  end do;
  eoms[general]:=tmp[theory_constrained];
  eprint(eoms,general,"So here is the general PGT once we removed all the zeroed params");
   #	a procedure to convert to observable quantities
  convert_variables:=proc(eq::equation,relations::set,old::function:=NULL)::equation;
    local unfinished,current,previous;
    current:=eq;
    unfinished:=true;
    while unfinished do
     previous:=current;
     current:=expand(subs(op(relations),current));
     unfinished:=evalb(simplify(previous)<>simplify(current));
    end do;
    simplify(current);
    if evalb(old=NULL) then
      return simplify(current);
    else
      return simplify(content(lhs(%),old)=0);
    end if;
  end proc;
  relations:=diff(a(t),t)=a(t)*H(t),Y(t)=a(t)*Q(t),X(t)=a(t)*U(t)/3-a(t)*H(t),Rho_d=sqrt(K)*rho_d(t)*a(t)^3,Rho_r=rho_r(t)*a(t)^4,Ls=K*rho_Ls,k=rho_k(t)*a(t)^2;
  for ii in X,Y,a,b do
    tmp[converted][ii]:=convert_variables(eoms[general][ii],{relations},a(t));
  end do;
  eprint(tmp,converted,"Initial attempt to convert to observable quantities");
  #	now we would like to introduce the omegas
  relations:=rho_r(t)=3*H(t)^2*O_r(t)/K,rho_d(t)=3*H(t)^2*O_d(t)/K,rho_Ls=3*H(t)^2*O_L(t)/K,rho_k(t)=-H(t)^2*O_k(t);
  for ii in X,Y,a,b do
    tmp[omegas][ii]:=convert_variables(tmp[converted][ii],{relations});
  end do;
  eprint(tmp,omegas,"Now we have put everything in the usual omega form");
  #	now for the shuffling 
  shuffle:=proc(form,eqs::table,version::name)::equation;
    local ii,wanted_coeffs,wanted_functions,wanted_coeffs_v,wanted_functions_v,labels,coef,buildvec,sys,sol,filler,n;
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
    filler:=[];
    n:=0;
    for ii in op(labels) do
      n:=n+1;
      filler:=[op(filler),ii=sol[n]];
    end do;
    return table([op(filler)]);
  end proc;
  form_H:=O_d(t)+O_r(t)+O_L(t);
  H_coeffs:=shuffle(form_H,tmp,omegas);
  form_q:=O_d(t)/2+O_r(t)-O_L(t);
  q_coeffs:=shuffle(form_q,tmp,omegas);
  initial_shuffle:=proc(form,resl)::equation;
    return simplify(solve(algsubs(form=x,simplify(resl[a]*tmp[omegas][a]+resl[b]*tmp[omegas][b])),x))=form;
  end proc;
  tmp[shuffle][H]:=initial_shuffle(form_H,H_coeffs);
  tmp[shuffle][q]:=initial_shuffle(form_q,q_coeffs);
  tmp[shuffle][X]:=tmp[omegas][X];
  tmp[shuffle][Y]:=tmp[omegas][Y];
  eprint(tmp,shuffle,"This is the very basic attempt at shuffling");
  eoms[Friedmann]:=tmp[shuffle];
  return NULL;
end proc;

tools[refine_eWGT_equations_of_motion]:=proc()::NULL;
  description "processing algorithm  for the eWGT equations";
  local tss,ii,tmp,composite_gauge_transform,specific_transform,convert_variables,relations,shuffle,form_H,H_coeffs,form_q,q_coeffs,initial_shuffle,simplifiers;
  global eoms;
  comment({bred,fyellow},"Now refining eWGT equations");
  #eoms:=table([]);
  tmp:=table([]);
  #	so first we get the raw versions
  eoms[raw][X]:=equations_of_motion[raw][1];
  eoms[raw][Y]:=equations_of_motion[raw][2];
  eoms[raw][a]:=equations_of_motion[raw][3];
  eoms[raw][Phi]:=equations_of_motion[raw][4];
  eprint(eoms,raw,"The raw equations of motion from the Lagrangian in terms of X, Y, a, Phi");
  #	next we want to get rid of any numerical factors
  for ii in X,Y,a,Phi do
    tmp[remove_prefactors][ii]:=simplify(eoms[raw][ii]/content(numer(normal(lhs(eoms[raw][ii])))));
  end do;
  #	we would also like to convert things into sigma notation 
  for ii in X,Y,a,Phi do
(*
    algsubs(6*A1+2*A2+2*A3+A4+A6=S3,tmp[remove_prefactors][ii],exact);
    algsubs(6*A1+A2+A3+A5-A6=S1,%,exact);
    algsubs(6*A1+2*A2+2*A3+3*A4-4*A5+A6=S2,%,exact);
    algsubs(2*A5-A4=(S3-S2)/2,%,exact);
    algsubs(B1+3*B2=U1,%,exact);
    algsubs(N=-6*U2,%,exact);
*)
    algsubs(6*A1+2*A2+2*A3+A4+A6=4*S3,tmp[remove_prefactors][ii],exact);
    algsubs(6*A1+A2+A3+A5-A6=4*S1,%,exact);
    algsubs(6*A1+2*A2+2*A3+3*A4-4*A5+A6=4*S2,%,exact);
    algsubs(2*A5-A4=2*(S3-S2),%,exact);
    algsubs(B1+3*B2=U1,%,exact);
    algsubs(N=-6*U2,%,exact);
    tmp[cosmological_parameters][ii]:=%
  end do;
  eprint(tmp,cosmological_parameters,"The equations of motion in a entirely general gauge with time, t, having removed any numerical prefactors and converted to sigma notation");
  #	now for the paper, we want to show that the cosmic theory parameters are well motivated
  for ii in X,Y,a,Phi do
    tss:=tmp[cosmological_parameters][ii];
    tss:=subs(b(t)=1,Phi(t)=1/sqrt(K),tss);
    tss:=simplify(subs(a(t)=a(t)*1,X(t)=Xw(t),Y(t)=Yw(t),%));
    tss:=dchange(t=f(u),tss,{u},known=f);
    tss:=evala(simplify(subs(diff(f(u),u)=a(u),tss)));
    tss:=simplify(subs(Xw(u)=X(u),Yw(u)=Y(u),tss));
    eoms[eWGT_conformal_form][ii]:=simplify(numer(normal(lhs(tss))))=0;
  end do;
  eprint(eoms,eWGT_conformal_form,"The conformal version of everything");
  simplifiers:={A,S1,S2,S3,U1,U2,N,L,Rho_r,Rho_d,k} minus {theory[Ricci]} minus {theory[quadratic_Riemann]} minus {theory[quadratic_torsion]} minus {theory[compensator_kinetic]} minus {theory[compensator_interaction]} minus {theory[radiation]} minus {theory[dust]} minus {theory[curvature]};
  for ii in X,Y,a,Phi do
    simplify(tmp[cosmological_parameters][ii],simplifiers);
    tmp[theory_constrained][ii]:=%;
  end do;
  eoms[general]:=tmp[theory_constrained];
  #	so this is a function which converts an equation of motion
  composite_gauge_transform:=proc(eq::equation)::equation;
    local tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    tmp1:=eq;
    tmp2:=subs(t=tp,tmp1);
    tmp3:=simplify(subs(X(tp)=Xpt(tp),Y(tp)=Ypt(tp),a(tp)=dtpdt(tp)*apt(tp),Phi(tp)=Phipt(tp)/dtpdt(tp),tmp2));
    tmp4:=simplify(dchange(tp=f(t),tmp3,{t},known=f));
    tmp5:=simplify(subs(diff(f(t),t)=dtpdt(t),tmp4));
    tmp6:=simplify(subs(Xpt(t)=Xp(t),Ypt(t)=Yp(t),apt(t)=ap(t),Phipt(t)=Phip(t),tmp5));
    return tmp6;
  end proc;
  for ii in X,Y,a,Phi do
    tmp[composite][ii]:=composite_gauge_transform(eoms[general][ii]);
  end do;
  eprint(tmp,composite,"Equations of motion in a general gauge, but expressed in terms of variables that come out of a generic composite transformation, the function dtpdt is dt'/dt, whilst t is the original time variable to be re-interpreted and all other fields are expressed with primes suppressed...");
  #	let us define a function which then uses the generically transformed equations and returns specific transformations, given the desired behavior of the fields in that gauge
  specific_transform:=proc(eq::name,required_behaviour::equation,timelabel::name)::equation;
    simplify(subs(Xp(t)=X(t),Yp(t)=Y(t),ap(t)=a(t),Phip(t)=Phi(t),tmp[composite][eq]));
    simplify(subs(required_behaviour,%));
    simplify(subs(t=timelabel,%));
    simplify(numer(normal(lhs(%)))=0);
    return %;
  end proc;  
  #	firstly the cosmic gauge 
  for ii in X,Y,a,Phi do
    tmp[cosmic][ii]:=specific_transform(ii,X(t)=-diff(a(t),t),t);
  end do;
  eoms[cosmic]:=tmp[cosmic];
  eprint(eoms,cosmic,"Equations of motion expressed in the variables Y,a,Phi in the cosmic gauge, with cosmic time, t");
  #	now for the conformal gauge
  for ii in X,Y,a,Phi do
    tmp[conformal][ii]:=specific_transform(ii,a(t)=a_0,x); 
  end do;
  #	and ridding ourselves of a_0
  for ii in X,Y,a,Phi do
    simplify(subs(N=N_s/a_0^2,L=L_s/a_0^4,Rho_d=Rho_d_s/a_0,tmp[conformal][ii])); 
    simplify(dchange(x=f(x_s),%,{x_s},known=f));
    simplify(subs(diff(f(x_s),x_s)=a_0,%));
    simplify(numer(normal(lhs(%)))=0);
    tmp[conformal][ii]:=%;
  end do;
  eoms[conformal]:=tmp[conformal];
  eprint(eoms,conformal,"Equations of motion expressed in the variables X,Y,Phi and a=a_0 in the conformal gauge, with dimensionless conformal time, x=eta/a_0, where eta is the conventional conformal time, and the dimensionful theory parameters N_s=N*a_0^2, L_s=L*a_0^4, Rho_d_s=Rho_d*a_0");
  #	now for the Einstein gauge
  for ii in X,Y,a,Phi do
    tmp[Einstein][ii]:=specific_transform(ii,Phi(t)=Phi_0,z); 
  end do;
  #	and ridding ourselves of Phi_0
  for ii in X,Y,a,Phi do
    simplify(subs(N=N_s/Phi_0^2,L=L_s/Phi_0^4,Rho_d=Rho_d_s/Phi_0,tmp[Einstein][ii])); 
    simplify(numer(normal(lhs(%)))=0);
    tmp[Einstein][ii]:=%;
  end do;
  eoms[Einstein]:=tmp[Einstein];
  eprint(eoms,Einstein,"Equations of motion expressed in the variables X,Y,a in the Einstein gauge, with Einstein time, z, and the dimensionful theory parameters N_s=N*Phi_0^2, L_s=L*Phi_0^4, Rho_d_s=Rho_d*Phi_0");
   #	a procedure to convert to observable quantities within the Einstein gauge 
  convert_variables:=proc(eq::equation,relations::set,old::function:=NULL)::equation;
    local unfinished,current,previous;
    current:=eq;
    unfinished:=true;
    while unfinished do
     previous:=current;
     current:=expand(subs(op(relations),current));
     unfinished:=evalb(simplify(previous)<>simplify(current));
    end do;
    simplify(current);
    if evalb(old=NULL) then
      return simplify(current);
    else
      return simplify(content(lhs(%),old)=0);
    end if;
  end proc;
  relations:=diff(a(z),z)=a(z)*H(z),Y(z)=a(z)*Q(z),X(z)=a(z)*U(z)/3-a(z)*H(z),Rho_d_s=rho_d_s(z)*a(z)^3,Rho_r=rho_r(z)*a(z)^4,L_s=rho_L_s,k=rho_k(z)*a(z)^2;
  for ii in X,Y,a,Phi do
    tmp[converted][ii]:=convert_variables(eoms[Einstein][ii],{relations},a(z));
  end do;
  eprint(tmp,converted,"Initial attempt to convert to observable quantities");
  #	now we would like to introduce the omegas
  relations:=rho_r(z)=3*H(z)^2*O_r(z)/K,rho_d_s(z)=3*H(z)^2*O_d(z)/K,rho_L_s=3*H(z)^2*O_L(z)/K,rho_k(z)=-H(z)^2*O_k(z);
  for ii in X,Y,a,Phi do
    tmp[omegas][ii]:=convert_variables(tmp[converted][ii],{relations});
  end do;
  eprint(tmp,omegas,"Now we have put everything in the usual omega form");
  return NULL;
end proc;

tools[walkthrough]:=proc(approach::name)::NULL;
  comment({byellow,underline,bold},"commencing walkthrough...");
  if evalb(approach=analytic) then
    if case=PGT_case_3 then
      read(`walkthroughs/cosh.mpl`);
    end if;
  elif evalb(approach=numerical) then
      read(`walkthroughs/D_numerical.mpl`);
  else
    error "you did not specify the walkthrough as being analytic or numerical, received %1", approach;
  end if;
  return NULL;
end proc;
