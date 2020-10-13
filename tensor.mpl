with(StringTools):
with(Physics):
Setup(mathematicalnotation=true):

rc:=proc()::NULL;
  description "resets the counters ii -- pp";
  global ii,jj,kk,ll,mm,nn,oo,pp;
  ii:='ii':
  jj:='jj':
  kk:='kk':
  ll:='ll':
  mm:='mm':
  nn:='nn':
  oo:='oo':
  pp:='pp':
  return NULL;
end proc;

index_translate:=proc(num::integer)
  description "translates numerical indices to the names of the spherical coordinates";
  if num=0 then
    return t;
  elif num=1 then
    return r;
  elif num=2 then
    return theta;
  elif num=3 then
    return phi;
  else
    error "expected numerical indices in the range 0--3";
  end if;
end proc;

sm:=proc(term,overind::set)
  description "performs the summation convention over a single term";
  local instances,tmp,ii,kk,jj;
  instances:={convert(term,string)};
  for ii in overind do
    tmp:={};
    for kk in instances do
      for jj from 0 to 3 do
	tmp:=tmp union {SubstituteAll(kk,convert(ii,string),convert(jj,string))};
      end do;
    end do;
    instances:=tmp;
  end do;
  tmp:=0;
  for ii in instances do
    tmp:=tmp+parse(ii);
  end do;
  tmp:=simplify(eval(tmp));
  return tmp;
end proc;

u:=proc(ii::integer)
  description "when we feed a counter or integer to a Physics vector, this procedure contravariantises the number with a tilde symbol";
  return parse(cat("~",convert(ii,string)));
end proc;
  
#	define the Lorentzian metric and Kroneker delta

for ii from 0 to 3 do
  for jj from 0 to 3 do
      etad||ii||d||jj:=g_[ii,jj];
      etau||ii||u||jj:=g_[u(ii),u(jj)];
      deltad||ii||u||jj:=g_[ii,u(jj)];
      deltau||ii||d||jj:=g_[u(ii),jj];
  end do;
end do;

rc():

for ii from 0 to 3 do
  for jj from 0 to 3 do
    for kk from 0 to 3 do
      for ll from 0 to 3 do
	  epsilonu||ii||u||jj||u||kk||u||ll:=LeviCivita[u(ii),u(jj),u(kk),u(ll)];
      end do;
    end do;
  end do;
end do;

rc();

#	definitions of the position vectors

Xu0:=t:
Xu1:=r*sin(theta)*cos(phi):
Xu2:=r*sin(theta)*sin(phi):
Xu3:=r*cos(theta):

#	the polar vectors

for jj in {t,r,theta,phi} do
  for ii from 0 to 3 do
    ed||jj||u||ii:=diff(Xu||ii,jj);
  end do;
end do;

rc();

#	the polar unit vectors

for jj in {t,r,theta,phi} do
  tmp:=simplify(sqrt(abs(sm(ed||jj||uii*ed||jj||ukk*etadiidkk,{ii,kk}))),trig) assuming r>0;
  for ll from 0 to 3 do
    ehd||jj||u||ll:=ed||jj||u||ll/tmp;
#    print(%);
  end do;
end do;

rc();

#	the components of the spin connection

for ii from 0 to 3 do
  for jj from 0 to 3 do
    #	the radial spin connection
    Au||ii||u||jj||dr:=sm(X*edtukk*edrull*(deltau||ii||dll*deltau||jj||dkk-deltau||ii||dkk*deltau||jj||dll)+Y*edtukk*edrull*sm(epsilonu||ii||u||jj||unnuoo*etadnndkk*etadoodll,{nn,oo}),{ll,kk})/sqrt(1-k*r^2);
    #	the polar spin connection
    Au||ii||u||jj||dtheta:=sm(edthetaukk*(X*edtull-(1/r)*(1-sqrt(1-k*r^2))*edrull)*(deltau||ii||dkk*deltau||jj||dll-deltau||jj||dkk*deltau||ii||dll)+Y*edtukk*edthetaull*sm(epsilonu||ii||u||jj||unnuoo*etadnndkk*etadoodll,{nn,oo}),{kk,ll});
    #	the azimuthal spin connection
    Au||ii||u||jj||dphi:=sm(edphiukk*(X*edtull-(1/r)*(1-sqrt(1-k*r^2))*edrull)*(deltau||ii||dkk*deltau||jj||dll-deltau||jj||dkk*deltau||ii||dll)+Y*edtukk*edphiull*sm(epsilonu||ii||u||jj||unnuoo*etadnndkk*etadoodll,{nn,oo}),{kk,ll});
  end do;
end do;

rc();

#	the PGT torsion

for ii from 0 to 3 do
  for jj from 0 to 3 do
    for kk from 0 to 3 do
      calTu||ii||d||jj||d||kk:=sm((-1/(b(t)*a(t)))*edtull*((X+diff(a(t)*b(t),t)/b(t))*(deltau||ii||d||jj*etadlld||kk-deltau||ii||d||kk*etadlld||jj)-(Y)*sm(epsilonu||ii||unnuooupp*etadnndll*etadood||jj*etadppd||kk,{nn,oo,pp})),{ll});
    end do;
  end do;
end do;

rc();

#	the PGT torsion invariants

beta_1:=sm(calTuiidjjdkk*calTulldnndmm*etadiidll*etaujjunn*etaukkumm,{ii,jj,kk,ll,nn,mm});
beta_2:=sm(calTuiidjjdkk*calTulldnndmm*deltadiiunn*deltaujjdll*etaukkumm,{ii,jj,kk,ll,nn,mm});
beta_3:=sm(calTuiidjjdkk*calTulldnndmm*deltadiiukk*deltadllumm*etaujjunn,{ii,jj,kk,ll,nn,mm});

