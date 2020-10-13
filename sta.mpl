#	This will form the basis of the sta package, we will write new procedures which replace those of the dirac package and work alongside

read `dirac/dirac_maple_test_r6.mpl`:
read `dirac/spherical_gravity_for_r10.mpl`:
with(dirac):

#	this just defines the Lorentz basis
gammad0:=gam0;
gammad1:=gam1;
gammad2:=gam2;
gammad3:=gam3;
gammau0:=gam0;
gammau1:=-gam1;
gammau2:=-gam2;
gammau3:=-gam3;

`type/multivector`:=proc(M)::bool;
  description "identifies multivectors as lists of 16 elements";
  evalb(type(M,list)=true and numelems(M)=16);
end proc;

`type/scalar`:=proc(M)::bool;
  description "identifies effective scalars";
  dirac[isscalar](M);
end proc;

`type/parsable`:=proc(x)::bool;
  description "identifies parsables as multivectors or effective scalars";
  evalb(dirac[ismulti](x)=true or dirac[isscalar](x)=true);
end proc;

sta[make_parsable]:=proc(x)::parsable;
  description "convert unparsable expressions to parsable ones, if possible, for use in products";
  if dirac[istricky](x) then
    ds(x)
  else
    x
  fi;
end proc;

sta[grade]:=proc(N,n::integer)::multivector;
  local M;
  description "projects out the required grade of a given multivector";
  M:=ds(N);
  if n<0 or n>4 then
    error "grade must be in the range 0-4 but received %1", n
  elif n=0 then
    ds(M[1]*scal)
  elif n=1 then
    dvpart(M)
  elif n=2 then
    dbpart(M)
  elif n=3 then
    dtpart(M)
  elif n=4 then
    ds(M[16]*ps)
  end if;
end proc;

#	here are the bootstrap versions of the extended products, requiring only the grade procedure.
(*
sta[`&!`]:=proc(x1,x2)::parsable;
  description "commutator product of two expressions";
  ds((x1&@x2-x2&@x1)/2);
end proc;

sta[`&.`]:=proc(x11,x22)::parsable;
  local tmp,x1,x2,ii,jj;
  description "interior product of two expressions";
  x1:=sta[make_parsable](x11);
  x2:=sta[make_parsable](x22);
  if type(x1,scalar) or type(x2,scalar) then
    x1&@x2
  else 
    tmp:=ds(0);
    for ii from 0 to 4 do
      for jj from 0 to 4 do
	tmp:=tmp+grade(grade(x1,ii)&@grade(x2,jj),abs(ii-jj));
      od;
    od;
    ds(tmp);
  end if;
end proc;

sta[`&^`]:=proc(x11,x22)::parsable;
  local tmp,x1,x2,ii,jj;
  description "exterior product of two expressions";
  x1:=sta[make_parsable](x11);
  x2:=sta[make_parsable](x22);
  if type(x1,scalar) or type(x2,scalar) then
    x1&@x2
  else 
    tmp:=ds(0);
    for ii from 0 to 4 do
      for jj from 0 to 4 do
	if ii+jj<=4 then
	  tmp:=tmp+grade(grade(x1,ii)&@grade(x2,jj),ii+jj);
	end if;
      od;
    od;
    ds(tmp);
  end if;
end proc;
*)

#	here are the fast extended products

sta[`&!`]:=proc(x11,x22)::parsable;
  local x1,x2;
  description "commutator product of two expressions";
  x1:=sta[make_parsable](x11);
  x2:=sta[make_parsable](x22);
  if type(x1,scalar) or type(x2,scalar) then
    0
  else 
    sta[commutator_product](x1,x2)
  end if;
end proc;

sta[`&.`]:=proc(x11,x22)::parsable;
  local tmp,x1,x2,ii,jj;
  description "interior product of two expressions";
  x1:=sta[make_parsable](x11);
  x2:=sta[make_parsable](x22);
  if type(x1,scalar) or type(x2,scalar) then
    x1&@x2
  else 
    sta[interior_product](x1,x2)
  end if;
end proc;

sta[`&^`]:=proc(x11,x22)::parsable;
  local tmp,x1,x2,ii,jj;
  description "exterior product of two expressions";
  x1:=sta[make_parsable](x11);
  x2:=sta[make_parsable](x22);
  if type(x1,scalar) or type(x2,scalar) then
    x1&@x2
  else 
    sta[exterior_product](x1,x2)
  end if;
end proc;

alias(grade=sta[grade]);

read `sta_products.mpl`:
