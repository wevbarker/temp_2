#	this is to better understand the Einstein freezing and parameter constraints

with(FileTools):
with(PDEtools):
with(CLIo);
with(ODEsys);

read(`Friedmann_form.mpl`):
tmp:=%:
for ii in X,Y,H,q do
  eoms[D_family][ii]:=tmp[ii]:
end do:

eprint(eoms,D_family,"These are the three equations for the whole of the cosmic class B family:"):

#	Now we would prefer to work with a system in which the U is eliminated

U_sol:=isolate(eoms[D_family][X],U(t)):

for ii in Y,H,q do
  tmp[simple][ii]:=simplify(subs(U_sol,eoms[D_family][ii])):
end do:

tmp[simple][auxiliary]:=tmp[simple][Y]:

#	Extract family of constant-Q solutions for radiation and dust

Q02_given_H:=proc(H_expr)::list;
  freezing_conditions:={diff(Q(t),t)=0,H(t)=H_expr}:
  frozen_auxiliary:=simplify(subs(freezing_conditions,tmp[simple][auxiliary])):
  frozen_einstein:=simplify(subs(freezing_conditions,tmp[simple][H])):
  frozen_auxiliary:=numer(normal(lhs(simplify(subs(Q(t)=Q0,K=1,frozen_auxiliary))))):
  frozen_einstein:=lhs(simplify(subs(Q(t)=Q0,K=1,frozen_einstein))):

  freezing_parameters:=proc(expr)
    local tmp:
    replace_S1:=simplify(subs(S1=S2*l,expr)):
    replace_S1_U3:=simplify(subs(U3=S2*w,replace_S1)):
    return %:
  end proc:

  frozen_auxiliary:=simplify(algsubs(Q0^2=Q02,simplify(freezing_parameters(frozen_auxiliary)/(16*S2^3*Q0*c)))):
  frozen_einstein:=simplify(algsubs(Q0^2=Q02,freezing_parameters(frozen_einstein))):

  Q02_solutions:=[solve(frozen_auxiliary,Q02)]:
  B_given_Q02:=proc(Q02v)
      return simplify(subs(w=U3/S2,simplify(subs(Q02=Q02v,frozen_einstein)))):
  end proc:
  B_solutions:=map(B_given_Q02,Q02_solutions):
  return [Q02_solutions,B_solutions]:
end proc:

radiation:=Q02_given_H(1/(2*t));
dust:=Q02_given_H(2/(3*t));
dark_energy:=Q02_given_H(HL);


B_d:=dust[2][2]:
B_L:=dark_energy[2][1]:
Q02_d:=dust[1][2]:
Q02_L:=dark_energy[1][1]:
simplify(B_d-B_L);
#this tells us 2-->1 constant-Q solution between matter and dark energy actually gives us the same B-factor!?

dit({},"and now for the radiation-dm solution");
B_r:=radiation[2][1]:
B_d:=dust[2][2]:
Q02_r:=radiation[1][1]:
Q02_d:=dust[1][2]:
crit:=simplify(numer(normal(subs(U3=1,B_d-B_r))));
solve(crit,l);
#	and so this tells us l=1 and l=1/2 are the only solutions for l here... I seem to remember both of these were pretty bad!

dit({},"Einstein modifiers");
subs(l=1.001,B_r):
print(%);
subs(l=1.001,B_d):
print(%);
subs(l=1.001,B_L):
print(%);
dit({},"Square torsion");
subs(l=1.001,Q02_r):
print(%);
subs(l=1.001,Q02_d):
print(%);
subs(l=1.001,Q02_L):
print(%);
#	okay, so there is a division by zero error at precisely l=1, but as we get close we can see that we get extremely close to the same B-factor in all cases, so long as we use negative U3.

(*	#this was used to show that there are no exact solutions.
pairtest:=proc(ii,jj)::NULL;
  local Q02_d,Q02_L:
  Q02_d:=simplify(subs(w=1,dust[1][ii])):
  Q02_L:=simplify(subs(w=1,dark_energy[1][jj])):
  B_d:=simplify(subs(U3=1,dust[2][ii])):
  B_L:=simplify(subs(U3=1,dark_energy[2][jj])):
  l_sols:=[solve(numer(B_d)=numer(B_L),l)]:
  print(%);
  print(ii,jj);
  for l_sol in l_sols do
    if l_sol=7/4 then
    Q02_d:=simplify(subs(l=l_sol,dust[1][ii])):
    Q02_L:=simplify(subs(l=l_sol,dark_energy[1][jj])):
    B_d:=simplify(subs(l=l_sol,dust[2][ii])):
    B_L:=simplify(subs(l=l_sol,dark_energy[2][jj])):
    print([l_sol,[Q02_d,Q02_L],[B_d,B_L]]);
  end if:
  end do:
  return NULL:
end proc:

for ii from 1 to 2 do
   for jj from 1 to 2 do
      pairtest(ii,jj):
  end do:
end do:

dit({},"last ditch for S2=1");
U3:=-1:
l:=-2:
w:=U3:
radiation;
dust;
dark_energy;
map(convert,radiation,float);
map(convert,dust,float);
map(convert,dark_energy,float);
*)








fin();
for ii from 1 to numelems(Q02_solutions) do
  Q02_r_||ii:=simplify(Q02_solutions[ii]);
  B_r_||ii:=simplify(subs(w=U3/S2,simplify(subs(Q02=Q02_r_||ii,frozen_einstein))));
end do;

dit({},"now try brute force for l which gives good B");

fin();













fin();
all_solutions:=[solve(frozen_auxiliary,Q0)];

fin();

find_Q02:=proc(fluid,instance::set)::set;
  local all_solutions,good_solutions,H_form,freezing_conditions,frozen_auxiliary;
  global l_radiation_1,l_radiation_2,l_radiation_3,l_radiation_4,l_radiation_5,l_matter_1,l_matter_2,l_matter_3,l_matter_4,l_matter_5,l_dark_energy_1,l_dark_energy_2,l_dark_energy_3,l_dark_energy_4,l_dark_energy_5;
  H_form:=table([radiation=1/(2*t),matter=2/(3*t),dark_energy=1,generic=c/t]);
  freezing_conditions:={diff(Q(t),t)=0,H(t)=eval(H_form[fluid])};
  frozen_auxiliary:=simplify(subs(freezing_conditions union instance,tmp[simple][auxiliary]));
  frozen_einstein:=simplify(subs(freezing_conditions union instance,tmp[simple][H]));
  frozen_auxiliary:=numer(normal(lhs(simplify(subs(Q(t)=Q0,K=1,frozen_auxiliary)))));
  frozen_einstein:=lhs(simplify(subs(Q(t)=Q0,K=1,frozen_einstein)));
  all_solutions:=[solve(frozen_auxiliary,Q0)];
  good_solutions:={}; 
  if instance={} then
    print(all_solutions);
    for ii from 1 to numelems(all_solutions) do
      l_||fluid||_||ii:=simplify(subs(Q0=all_solutions[ii],frozen_einstein)); 
    end do;
    Q_solutions_||fluid:=all_solutions;
    good_solutions:=%;
  else
    for ii in all_solutions do
      trial:=convert(ii,float)^2;
      if trial>0 then
	good_solutions:=good_solutions union {log(abs(convert(ii,float)))};
      end if;
    end do;
  end if;
  return good_solutions;
end proc:

find_Q02(generic,{});
tmp:=%;
sol:=tmp[2];
dit({},"here is the generaic solution for c/t=H");
print(sol);
replace_S1:=simplify(subs(S1=S2*l,sol));
replace_U3:=simplify(subs(U3=S2*w,replace_S1));
simplify(%,radical);
simplify(%^2) assuming real;

fin();

find_Q02(radiation,{});
print(%);
find_Q02(matter,{});
print(%);
find_Q02(dark_energy,{});
print(%);

dit({},"matter to dark energy smooth transition");

for ii from 1 to numelems(Q_solutions_matter) do
  for jj from 1 to numelems(Q_solutions_dark_energy) do
     print(ii,jj);
     solutions:=solve({l_matter_||ii=l_dark_energy_||jj},{S1});
     for kk in solutions do
	dit({},"solution...");
	print(kk);
	print(simplify(subs(kk,S2=1,U3=1,l_matter_||ii)));
     end do;
  end do;
end do;

dit({},"radiation to matter smooth transition");

for ii from 1 to numelems(Q_solutions_radiation) do
  for jj from 1 to numelems(Q_solutions_matter) do
     print(ii,jj);
     solutions:=solve({l_radiation_||ii=l_matter_||jj},{S1});
     for kk in solutions do
	dit({},"solution");
	print(kk);
	print(simplify(subs(kk,S2=1,U3=1,l_matter_||jj)));
     end do;
  end do;
end do;


fin();
