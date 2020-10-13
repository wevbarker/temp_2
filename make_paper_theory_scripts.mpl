read(`theory_tools.mpl`):

grand:=table([]);

for ii from 1 to 33 do
  grand[case||ii]:=0;
end do;

for ii from 1 to 33 do
  tmp:=testcase(ii):
  caselist:=op(map(op,{indices(grand)})):
  original:=true:
  for jj in caselist do
    if evalb(tmp[3]=grand[jj][constraints]) then
      original:=false;
      caseis:=jj;
    end if;
  end do;
  if original then
    nm:=case||ii;
    grand[nm][constraints]:=tmp[3];
    if tmp[2]=1 then
      grand[nm][theories]:={tmp[1]};
      grand[nm][theories2]:={};
      grand[nm][theories3]:={};
    elif tmp[2]=2 then
      grand[nm][theories]:={};
      grand[nm][theories2]:={tmp[1]};
      grand[nm][theories3]:={};
    elif tmp[2]=3 then
      grand[nm][theories]:={};
      grand[nm][theories2]:={};
      grand[nm][theories3]:={tmp[1]};
    end if;
  else
    if tmp[2]=1 then
      grand[caseis][theories]:=grand[caseis][theories] union {tmp[1]};
    elif tmp[2]=2 then
      grand[caseis][theories2]:=grand[caseis][theories2] union {tmp[1]};
    elif tmp[2]=3 then
      grand[caseis][theories3]:=grand[caseis][theories3] union {tmp[1]};
    end if;
  end if;
end do:
dit({fblack,byellow,bold,underline},"Here is the raw cosmic theory table:");
print(grand);
#W();
great:=table([]):

for jj from 1 to 10 do
  caseindex:=1;
  for ii from 1 to 33 do
    if evalb(grand[case||ii]<>0) then
      if evalb(numelems(grand[case||ii][constraints])=jj) then
	great[O||jj||N||caseindex]:=eval(grand[case||ii]);
	caseindex:=caseindex+1;
      end if;
    end if;
  end do;
end do:

dit({},"Case table ordered by number of cosmic constraints:");
print(great);

#(*
great[O3N||6][constraints]:={S2-S3,U2,cA};
great[O3N||6][theories]:={};
great[O3N||6][theories2]:={};
great[O3N||6][theories3]:={};
#*)

k:=1;
for jj from 1 to 10 do
  for ii from 1 to 10 do
    if ii<>jj and assigned(great[O||3||N||jj]) and assigned(great[O||3||N||ii]) then
      inter:=great[O||3||N||jj][constraints] intersect great[O||3||N||ii][constraints];
      if evalb(numelems(inter)=2) then
	original:=true;
	for kk from 1 to 10 do
	  if great[O||2||N||kk][constraints]=inter then
	    original:=false;
	  end if;
	end do;
	if original then
	  great[O2N||k][constraints]:=inter;
	  great[O2N||k][theories]:={};
	  great[O2N||k][theories2]:={};
	  great[O2N||k][theories3]:={};
	  k:=k+1;
	end if;
      end if;
    end if;
  end do;
end do;    

k:=1;
for jj from 1 to 10 do
  for ii from 1 to 10 do
    if ii<>jj and assigned(great[O||2||N||jj]) and assigned(great[O||2||N||ii]) then
      inter:=great[O||2||N||jj][constraints] intersect great[O||2||N||ii][constraints];
      if evalb(numelems(inter)=1) then
	original:=true;
	for kk from 1 to 10 do
	  if great[O||1||N||kk][constraints]=inter then
	    original:=false;
	  end if;
	end do;
	if original then
	  great[O1N||k][constraints]:=inter;
	  great[O1N||k][theories]:={};
	  great[O1N||k][theories2]:={};
	  great[O1N||k][theories3]:={};
	  k:=k+1;
	end if;
      end if;
    end if;
  end do;
end do;    



paternity_test:=proc(parent::name,child::name)
  description "determines if one cosmic class is a child of another via a single constraint, returns [true/false,{single constraint}]";
  global great;
  if verify(great[parent][constraints],great[child][constraints],`subset`) then
    return [true,great[child][constraints] minus great[parent][constraints]];
  elif evalb(numelems(solve(simplify(great[child][constraints],great[parent][constraints])))=1) then
    return [true,{lhs(op(solve(simplify(great[child][constraints],great[parent][constraints]))))}];
  else
    dit({},"allegedly really hard for maple...");
    par:=great[parent][constraints]:
    chi:=great[child][constraints]:
    commonality:=par intersect chi:
    par:=par minus commonality:
    chi:=chi minus commonality:
    dit({},"stripped-down versions...");
    if numelems(map(primpart,simplify(chi,solve(par))))=1 then
      return [true,{op(map(primpart,simplify(chi,solve(par))))}]:
    else
      return [false,{}];
    end if:
  end if;
end proc;
  
for ii from 1 to 5 do
  for jj from 1 to 10 do
    child:=O||ii||N||jj;
    if assigned(great[child]) then
      par:={};
      for kk from 1 to 10 do
	parent:=O||(ii-1)||N||kk;
	if assigned(great[parent]) then
	  if paternity_test(parent,child)[1] then
	    par:=par union {parent};
	  end if;
	end if;
      end do;
      great[child][parents]:=par;
    end if;
  end do;
end do:

for ii in indices(great) do
  if great[ii][parents]={} then
    print(great[ii]);
  end if;
end do;

great[O0N1][constraints]:={}:
great[O0N1][parents]:={}:
great[O0N1][theories]:={}:
great[O1N1][parents]:={O0N1}:

for ii in indices(great) do
  great[ii[1]][alltheories]:=great[ii[1]][theories] union great[ii[1]][theories2] union great[ii[1]][theories3]:
end do:

checkifin:=proc(node::name,theory::integer)::boolean;
  global great;
  tmmp:=eval(great[node][alltheories]):
  if type(tmmp,set) then
    if paper_to_intermediate[theory] in tmmp then
      return true;
    else
      return false;
    end if:
  else
    return false;
  end if:
end proc:

print(great[O0N1]);

isstrong:=proc(parent::name,child::name)::boolean;
  global great;
  if checkifin(parent,2) and checkifin(child,1) then
    return true;
  elif checkifin(parent,2) and checkifin(child,15) then
    return true;
  elif checkifin(parent,2) and checkifin(child,16) then
    return true;
  elif checkifin(parent,15) and checkifin(child,14) then
    return true;
  elif checkifin(parent,15) and checkifin(child,12) then
    return true;
  elif checkifin(parent,16) and checkifin(child,14) then
    return true;
  elif checkifin(parent,16) and checkifin(child,11) then
    return true;
  elif checkifin(parent,1) and checkifin(child,12) then
    return true;
  elif checkifin(parent,1) and checkifin(child,11) then
    return true;
  elif checkifin(parent,14) and checkifin(child,10) then
    return true;
  elif checkifin(parent,12) and checkifin(child,10) then
    return true;
  elif checkifin(parent,11) and checkifin(child,10) then
    return true;
  else
    return false;
  end if:
end proc:

#fin();


letters:=table([1=A,2=B,3='C',4='D',5='E',6=F,7=G,8=H,9='I',10=J,11=K,12=L,13=N,14=M,15=O,16=P,17=Q,18=R]);
#	in order to get the paper version, we use black, and for the poster we use purple...
cosmic_coefficients:=table([cA="{\\color{black}\\check{\\alpha}_0}",S1="{\\color{black}\\sigma_1}",S2="{\\color{black}\\sigma_2}",S3="{\\color{black}\\sigma_3}",U1="{\\color{black}\\upsilon_1}",U2="{\\color{black}\\upsilon_2}"]);

name_generator:=proc(expr::name,num::integer)::string;
  global count, cosmic_class_names;
  print(count);
  level:=numelems(great[expr][constraints]);
  if evalb(expr=O0N1) then
    nameis:="PGT/eWGT";
  else
    cons:=great[expr][constraints];
    YM:=false;
    ks:=false;
    m:=false;
    if evalb(cA in cons) then
      YM:=true;
    end if;
    if evalb(cA in cons) and evalb(S3 in cons) then
     ks:=true;
    end if;
    if evalb(U1 in cons)=false or evalb(U2 in cons)=false then
      m:=true;
    end if;
    sYM:=proc(a::boolean)::string;
      if a then
	#return "\\scriptstyle{q}";
	return "\\phantom{\\scriptstyle{q}}";
      else
	return "\\phantom{\\scriptstyle{q}}";
      end if;
    end proc;
    sks:=proc(a::boolean)::string;
      if a then
	return "\\scriptstyle{\\slashed{k}}";
      else
	return "\\phantom{\\scriptstyle{\\slashed{k}}}";
      end if;
    end proc;
    sm:=proc(a::boolean)::string;
      if a then
	return "\\phantom{*}";
      else
	return "\\phantom{*}";
      end if;
    end proc;
    prop:=cat("$\\substack{",sm(m),"\\\\",sks(ks),"\\\\",sYM(YM),"}$");
    colr:=proc(bl::boolean)::string;
      if bl then
	return "{\\color{goodcosmic}";
      else
	return "{\\color{badcosmic}";
      end if;
    end proc;
    col:=colr(m);
    nameis:=cat("\\textsuperscript{",convert(level,string),"}",col,convert(letters[count],string),prop,"}");
    #cosmic_class_names[expr]:=convert(letters[count],string):
    cosmic_class_names[expr]:=cat("\\textsuperscript{",convert(level,string),"}",convert(letters[count],string)):
  end if;
  return nameis;
end proc;

label_generator:=proc(expr::name,num::integer)::string;
  level:=numelems(great[expr][constraints]);
  nameis:=cat(convert(level,string),convert(letters[num],string));
  return nameis;
end proc;

QFT_case_generator:=proc(cas::integer)::string;
  properties:=theory_properties[cas];
  if properties={Mv,Ml,s2p} then
    prop:="$\\substack{2^+\\\\\\circ\\circ\\\\\\bullet}$";
  elif properties={Mv,s2p} then
    prop:="$\\substack{\\phantom{2^+}\\\\\\vphantom{\\circ\\circ}\\\\\\bullet}$";
  elif properties={Ml,s2p} then
    prop:="$\\substack{2^+\\\\\\circ\\circ\\\\\\vphantom{\\bullet}}$";
  elif properties={Mv} then
    prop:="$\\substack{\\vphantom{2^+}\\\\\\vphantom{\\circ\\circ}\\\\\\bullet}$";
  elif properties={Ml} then
    prop:="$\\substack{\\vphantom{2^+}\\\\\\circ\\circ\\\\\\vphantom{bullet}}$";
  end if;
  if numelems(eval(case_translation[cas]))=1 then
    if cas in {4,5,6,8} then
      return cat("{\\color{goodquantum}",convert(case_translation[cas][1],string),prop,"}");
    else
      return cat("{\\color{badquantum}",convert(case_translation[cas][1],string),prop,"}");
    end if:
  else
    prop1:=cat("$^{*",convert(case_translation[cas][2],string),"}$"):
    if cas in {4,5,6,8} then
      return cat("{\\color{goodquantum}",prop1,convert(case_translation[cas][1],string),prop,"}");
    else
      return cat("{\\color{badquantum}",prop1,convert(case_translation[cas][1],string),prop,"}");
    end if:
  end if:
end proc;

QFT_label_generator:=proc(theories::set)::string;
  theoriess:=sort(convert(theories,list));
  label:=QFT_case_generator(theoriess[1]);
  if numelems(theoriess)>1 then
    for ii from 2 to numelems(theoriess) do
      label:=cat(label,",",QFT_case_generator(theoriess[ii]));
    end do;
  end if;
  return label;
end proc;

print_cosmiccase:=proc(row::integer,x,y,lab::string,nam::string)::NULL;
  global count,cosmic_class_names;
  namm:=convert(lab,name);
  theoriess:=eval(great[namm][theories]);
  theoriess2:=eval(great[namm][theories2]);
  theoriess3:=eval(great[namm][theories3]);
  all:=theoriess union theoriess2 union theoriess3;
  isparent:=evalb(namm in parentlist[row]);
  X1:=x+8/3;
  Y1:=y-1/7;
  X2:=x+8/3;
  Y2:=y-1/2;
  X3:=x+8/3;
  Y3:=y-7/4;
  
  if row=0 then
    Yd:=y-2/3;
  elif row<=1 then
    Yc:=y+3/3;
    Yd:=y-3/2;
  elif row=2 then
    Yc:=y+3/3;
    Yd:=y-3/2;
  elif row>2 then
    Yc:=y+3/3;
    Yd:=y-3/2;
  end if;
  trow||row||upper:=Yc;
  trow||row||lower:=Yd;
  if evalb(lab="O0N1") then
      printf("\\node [cosmiccase] at (%f,%f) (%s) {%s};\n",x,y,lab,"PGT\\textsuperscript{+,q}");
      printf("\\node [tie] at (%f,%f) (%s) {};\n",x,Yd,cat(lab,"d"));
      printf("\\draw [fill=black] ($(%s.south)+(-1/16,0)$) -- (%s.center) -- ($(%s.south)+(1/16,0)$) -- cycle;\n",lab,cat(lab,"d"),lab);
  elif evalb(all={}) then
      printf("\\node [tie] at (%f,%f) (%s) {};\n",x,Yc,cat(lab,"c"));
      printf("\\node [cosmiccase] at (%f,%f) (%s) {%s};\n",x,y,lab,"$.$");
      printf("\\node [tie] at (%f,%f) (%s) {};\n",x,Yd,cat(lab,"d"));
      printf("\\draw [fill=black] (%s.center) -- ($(%s.center)+(-1/16,0)$) -- (%s.center) -- ($(%s.center)+(1/16,0)$) -- cycle;\n",cat(lab,"c"),lab,cat(lab,"d"),lab);
  else
      printf("\\node [tie] at (%f,%f) (%s) {};\n",x,Yc,cat(lab,"c"));
      printf("\\node [cosmiccase] at (%f,%f) (%s) {%s};\n",x,y,lab,nam);
      printf("\\node [tie] at (%f,%f) (%s) {};\n",x,Yd,cat(lab,"d"));
      printf("\\draw [fill=black] (%s.center) -- ($(%s.north)+(-1/16,0)$) -- ($(%s.north)+(1/16,0)$) -- cycle;\n",cat(lab,"c"),lab,lab);
      if isparent then
	printf("\\draw [fill=black] (%s.center) -- ($(%s.south)+(-1/16,0)$) -- ($(%s.south)+(1/16,0)$) -- cycle;\n",cat(lab,"d"),lab,lab);
      end if;
      count:=count+1:
      if evalb(theoriess={})=false then
	printf("\\node [QFTcase] at (%f,%f) (%s) {%s};\n",X1,Y1,cat(lab,"QFT1"),QFT_label_generator(theoriess));
	printf("\\path [rel,densely dotted] (%s) -- (%s.west);\n",lab,cat(lab,"QFT1"));
      end if;
      if evalb(theoriess2={})=false then
	printf("\\node [QFTcase] at (%f,%f) (%s) {%s};\n",X2,Y2,cat(lab,"QFT2"),QFT_label_generator(theoriess2));
	printf("\\path [rel,densely dotted,postaction={on each segment={PlusTwoNonCosmic}}] (%s) -- (%s.west);\n",lab,cat(lab,"QFT2"));
      end if;
      if evalb(theoriess3={})=false then
	printf("\\node [QFTcase] at (%f,%f) (%s) {%s};\n",X3,Y3,cat(lab,"QFT3"),QFT_label_generator(theoriess3));
	printf("\\path [rel,densely dotted,postaction={on each segment={PlusThreeNonCosmic}}] (%s) -- (%s.west);\n",lab,cat(lab,"QFT3"));
      end if;
  end if;
  return NULL;
end proc;

print_constraint:=proc(row::integer,parent::name,child::name,x,y)::NULL;
  ii:=substring(convert(parent,string),-1);
  middle:=1/200;
  if assigned(label_shift[parent][child]) then
    top:=eval(label_shift[parent][child]):
    middle:=eval(label_shift[parent][child]):
    bottom:=eval(label_shift[parent][child]):
  end if:
  constraint:=op(paternity_test(parent,child)[2]);
  if constraint<>0 then
    constraint:=sign(convert(constraint,list)[1])*constraint;
  end if;
  consr:=convert(constraint,string);
  for ii in cA,S1,S2,S3,U1,U2 do
    str:=eval(cosmic_coefficients[ii]);
    consr:=Substitute(consr,convert(ii,string),str);
  end do;
  star:=proc(s::string)::boolean;
    return evalb(s="*");
  end proc;
  consr:=Remove(star,consr);
  consr:=cat("${",consr,"}$");
  consname:=cat(convert(parent,string),convert(child,string));
    parnam:=convert(parent,string);
    chinam:=convert(child,string);
  #if type(constraint,name) then
    printf("\\path [name path=%s] (-20,%f) -- (15,%f);\n",cat(consname,"aid"),y+middle,y+middle);
(*
    print(parent);
    print(child);
    print(great[parent]);
    print(great[child]);
    print(isstrong(parent,child));
    fin();
*)
    if isstrong(parent,child) then
      printf("\\draw [line width=0.8mm,name path=%s] (%s.center) to[out=-90,in=90,loop,looseness=0.4]  (%s.center);\n",cat(parnam,chinam),cat(parnam,"d"),cat(chinam,"c"));
    else
      printf("\\draw [line width=0.4mm,name path=%s] (%s.center) to[out=-90,in=90,loop,looseness=0.4]  (%s.center);\n",cat(parnam,chinam),cat(parnam,"d"),cat(chinam,"c"));
    end if:
    printf("\\draw[name intersections={of=%s and %s}] node at (intersection-1) [fill=white] {%s};\n",cat(parnam,chinam),cat(consname,"aid"),consr);
  return NULL;
end proc;

rownums:=table([]);
constraintnums:=table([]);
for row from 0 to 5 do
  n:=0:
  for ii from 1 to 10 do
    if assigned(great[O||row||N||ii]) then
      n:=ii;
    end if;
  end do;
  rownums[row]:=n:
  m:=0;
  for kk from 1 to n do
   m:=m+numelems(great[O||row||N||kk][parents]);
  end do;
  constraintnums[row]:=m
end do;

parentlist:=table([5={}]);
redo_parentlist:=proc()::NULL;
  global parentlist;
  for ii from 0 to 4 do
    lis:={};
    row:=ii+1;
    for jj from 1 to rownums[row] do
      lis:=lis union great[O||row||N||jj][parents];
    end do;
    parentlist[ii]:=lis;
  end do;
end proc;

cosmic_x:=proc(row::integer,ii::integer)::float;
  local n,shift4,shift5;
  n:=rownums[row];
  if row=2 then
    return (3*42/7)*(ii-(3)/2-1/6); 
  elif row=1 or row=0 then
    return -(2*42/7); 
  elif row=5 then
    return (2*42/7)*(ii-(3+1)/2); 
  elif row=3 then
    return (42/7)*(ii-(5+1)/2)-42/7; 
  elif row=4 then
    return (42/7)*(ii-(5+1)/2)-42/7; 
  else
    return (42/n)*(ii-(n+1)/2); 
  end if;
end proc;

constraint_x:=proc(row::integer,k::integer)::float;
 m:=constraintnums[row];
  return (20/m)*(k-(m+1)/2);
end proc;

cosmic_y:=proc(row::integer)::float;
    if row<2 then
      y:=-2.5*row;
    elif row>=2 then
      y:=-2.5-8*(row-1);
    end if;
  return y;
end proc;

constraint_y:=proc(row::integer)::float;
  if row<=3 then
    y:=cosmic_y(row);
  elif row>3 then
    y:=cosmic_y(row)+2;
  end if;
  return y;
end proc;

rowprint:=proc(row::integer)::NULL;
  global count,cosmic_class_names;
  k:=1:
  for ii from 1 to rownums[row] do
    cas:=O||row||N||ii;
    print_cosmiccase(row,cosmic_x(row,ii),cosmic_y(row),convert(cas,string),name_generator(cas,ii));
    pars:=convert(great[cas][parents],list);
    if evalb(assigned(great[cas][parents])) then
      for jj in op(pars) do
	parnam:=convert(jj,string);
	urow:=row-1;
	print_constraint(row,jj,cas,constraint_x(row,k),(1/2)*(trow||row||upper+trow||urow||lower));
	k:=k+1;
      end do;
    end if;
  end do;
  return NULL;
end proc:

count:=0;
cosmic_class_names:=table([]);
recast:=proc()::NULL;
  global count,parentlist,cosmic_class_names;
  redo_parentlist();
  dit({},"recasting...");
  count:=1;
  for ii from 0 to 5 do
    rowprint(ii);
  end do;
  count:=1;
  writeto(`/home/williamb/Documents/physics/papers/paper-2/paper/aps/graph.tex`);
  for ii from 0 to 5 do
    rowprint(ii);
  end do;
  writeto(terminal);
  count:=1;
  for ii from 0 to 5 do
    rowprint(ii);
  end do;
end proc;

recast();

subsor:=proc(targ,rel)
 return subs(rel,targ) ;
end proc;

nodeswap:=proc(row::integer,ii::integer,jj::integer)::NULL;
    global great,cosmic_class_names;
    A:=O||row||N||ii;
    B:=O||row||N||jj;
    dit({},"swapping a node...");
    tmp:=eval(great[B]);
    great[B]:=eval(great[A]);
    great[A]:=eval(tmp);
    if row<>5 then
      for ll from row+1 to 5 do
	for kk from 1 to rownums[ll] do
	  great[O||ll||N||kk][parents]:=map(subsor,great[O||ll||N||kk][parents],A=tmp2);
	  great[O||ll||N||kk][parents]:=map(subsor,great[O||ll||N||kk][parents],B=A);
	  great[O||ll||N||kk][parents]:=map(subsor,great[O||ll||N||kk][parents],tmp2=B);
	end do;
      end do;
    end if;
    recast();
    return NULL;
end proc;

print(eval(great[O4N4]));

nodeswap(3,1,4);
nodeswap(3,1,3);
nodeswap(4,1,5);
nodeswap(3,3,4);
nodeswap(3,5,6);
nodeswap(3,4,5);
nodeswap(4,4,5);
nodeswap(3,1,2);
nodeswap(4,2,3);

nodeswap(4,1,2);
nodeswap(4,2,3);
nodeswap(4,3,4);
nodeswap(4,4,5);
nodeswap(4,1,3);
nodeswap(3,1,2);

label_shift:=table([]):
label_shift[O3N2][O4N4]:=-3/2:
label_shift[O3N5][O4N7]:=-3/2:

label_shift[O2N2][O3N3]:=1:
label_shift[O3N4][O4N5]:=1:

label_shift[O3N1][O4N2]:=-1:
label_shift[O3N2][O4N1]:=1:
label_shift[O3N2][O4N3]:=1:
label_shift[O3N3][O4N2]:=-1:

label_shift[O3N5][O4N5]:=1:
label_shift[O3N4][O4N4]:=-1:
recast();

#	this produces macros for the paper, for critical cases and cosmic classes

cosmic_class_list:=table([]):
for kk in indices(case_translation) do
  for jj in indices(great) do
    if assigned(great[jj[1]][theories]) and kk[1] in great[jj[1]][theories] then
      cosmic_class_list[kk[1]]:=eval(cosmic_class_names[jj[1]]):
    elif assigned(great[jj[1]][theories2]) and kk[1] in great[jj[1]][theories2] then
      cosmic_class_list[kk[1]]:=eval(cosmic_class_names[jj[1]]):
    elif assigned(great[jj[1]][theories3]) and kk[1] in great[jj[1]][theories3] then
      cosmic_class_list[kk[1]]:=eval(cosmic_class_names[jj[1]]):
    end if:
  end do: 
end do:
print(indices(great));
print(eval(cosmic_class_list));
print(eval(cosmic_class_names));
set_of_cosmic_class_names:={}:
for ii in indices(case_translation) do
  set_of_cosmic_class_names:=set_of_cosmic_class_names union {cosmic_class_list[ii[1]]}:
end do:
number_of_cosmic_classes:=numelems(set_of_cosmic_class_names):

writeto(`/home/williamb/Documents/physics/papers/paper-2/paper/aps/cosmic_class_names.tex`):
printf("\\newrobustcmd{\\numbercosmicclass}{%d}%%\n",number_of_cosmic_classes):
printf("\\newrobustcmd{\\cosmicclass}[1]{%%\n"):
printf("\\IfEqCase{#1}{%%\n"):
for ii in indices(case_translation) do
  printf("{%s}{Class %s}%%\n",convert(case_translation[ii[1]][1],string),cat(cosmic_class_list[ii[1]],"\\xspace")):
end do:
printf("{null}{Class %s}%%\n",cat(cosmic_class_list[11],"*\\xspace")):
printf("}[\\packageError{cosmicclass}{Unidentified Critical Case: #1}{}]%%\n"):
printf("}\n"):
printf("\\newrobustcmd{\\criticalcase}[1]{%%\n"):
printf("\\IfEqCase{#1}{%%\n"):
for ii in indices(case_translation) do
  if numelems(case_translation[ii[1]])=1 then
    printf("{%s}{Case %s}%%\n",convert(case_translation[ii[1]][1],string),cat(convert(case_translation[ii[1]][1],string),"\\xspace")):
  else
    printf("{%s}{Case %s}%%\n",convert(case_translation[ii[1]][1],string),cat("\\textsuperscript{*",convert(case_translation[ii[1]][2],string),"}",convert(case_translation[ii[1]][1],string),"\\xspace")):
  end if:
end do:
printf("}[\\packageError{criticalcase}{Unidentified Critical Case: #1}{}]%%\n"):
printf("}\n"):
writeto(terminal):

print(great);

fin();
