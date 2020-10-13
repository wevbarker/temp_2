with(StringTools):
with(LinearAlgebra):
with(CLIo):
read `tools.mpl`:
with(tools):

read(`parameters.mpl`):
read(`particle_content.mpl`):

(*
particle_content:=table([1=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),2=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),3=table([J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[Al2],J2Pm=[Al2]]),4=table([J1Pm=[Al2],J1Pp=[Al2],J2Pm=[Al2]]),5=table([J1Pm=[AL2],J1Pp=[Al2],J2Pp=[Al2]]),6=table([J0Pm=[Al2],J1Pm=[Al2],J1Pp=[Al2],J2Pp=[Al2]]),7=table([J0Pm=[Al0],J1Pm=[Al2],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),8=table([J0Pp=[Al2],J1Pm=[Al2],J1Pp=[Al2],J2Pm=[Al2]]),9=table([J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[Al2],J2Pp=[Al2]]),10=table([J0Pm=[Al0],J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),11=table([J0Pm=[Al2],J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[Al2],J2Pp=[Al2]]),12=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al0,Al0],[Al0,sl2],[Al0,al2]],J1Pp=[[Al0,Al0],[Al0,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),13=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al0,Al0],[Al0,sl2],[Al0,al2]],J1Pp=[[Ali,Al0],[Ali,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),14=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Ali,Al0],[Ali,sl2],[Ali,al2]],J1Pp=[[Al0,Al0],[Al0,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),15=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al0,Al0],[Al0,sl2],[Al0,al2]],J1Pp=[[Al0,Al0],[Al0,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),16=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[[Al2,Al0],[Al2,al2]]]),17=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[Al0,sl2,al2],J1Pp=[Al0,al2]]),18=table([J0Pm=[Av2],J1Pp=[Al0,al2]]),19=table([J0Pm=[Av2],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),20=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[[Al2,Al0],[Al2,al2]]]),21=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[Al0,al2],J2Pm=[Al2]]),22=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[Al0,al2],J2Pm=[Al2]]),23=table([J0Pm=[Av2],J0Pp=[Al2],J1Pp=[Al0,al2]]),24=table([J0Pm=[Av2],J1Pm=[Al2,sl2,al2],J1Pp=[[Al0,al2],[Al0,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),25=table([J0Pm=[Av2],J1Pm=[Al2,sl2,al2],J1Pp=[[Ali,al2],[Ali,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),26=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al2],J1Pp=[[Al2,Al0],[Al2,al2]]]),27=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[Al0,al2],J2Pp=[Al2]]),28=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al2],J1Pp=[Al0,al2],J2Pm=[Al2]]),29=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[[Al2,Al0],[Al2,sl2],[Al2,al2]],J1Pp=[Al0,al2],J2Pp=[Al2]]),30=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al0,sl2,al2],J1Pp=[[Al0,Al0],[Al0,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),31=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al0,sl2,al2],J1Pp=[[Ali,Al0],[Ali,al2]],J2Pm=[Al0],J2Pp=[Al0,sl2]]),32=table([J0Pm=[Av2],J0Pp=[Al2],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),33=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[Al0,sl2,al2],J1Pp=[[Al2,Al0],[Al2,al2]],J2Pp=[Al2]]),34=table([J0Pm=[Av2],J0Pp=[sv2],J1Pm=[av2],J1Pp=[Al2],J2Pm=[sl2],J2Pp=[al2]]),35=table([J0Pm=[[Av2,sv2]],J0Pp=[sv2,av2],J1Pm=[av2],J1Pp=[Al2],J2Pm=[sl2],J2Pp=[al2]])]);

parameter_constraints:=table([1={R1,R3/2-R4,T1,T3,L},2={R1,R3/2-R4,T1,L},3={R2,R1-R3,R4,T1,T2,L},4={R2,R1-R3,R4,T1,T2,T3,L},5={R1,R2,R3/2-R4,T1,T2,T3,L},6={R1,R3/2-R4,T1,T2,T3,L},7={R1,R2,R3/2-R4,T1,T3,L},8={R2,2*R1-2*R3+R4,T1,T2,T3,L},9={R1,R2,R3/2-R4,T1,T2,L},10={R1,R2,R3/2-R4,T1,L},11={R1,R3/2-R4,T1,T2,L},12={R1,R3,R4,R5,L},13={R1,R3,R4,R5,T1+T2,L},14={R1,R3,R4,T1+T3,L},15={R1,R3,R4,R5,T1+T2,T1+T3,L},16={R1,R3,R4,T1,L},17={R1,R3,R4,R5,T1,L},18={R1,R3,R4,R5,T1,T3,L},19={R1,R3/2-R4,R3/2-R5,T1,T3,L},20={R1,R3,R4,T1,T3,L},21={R1-R3,R4,2*R1+R5,T1,L},22={R1-R3,R4,2*R1+R5,T1,T3,L},23={R1,2*R3-R4,2*R3+R5,T1,T3,L},24={R1,R3,R4,R5,T3,L},25={R1,R3,R4,R5,T1+T2,T3,L},26={R1,2*R3-R4,T1,T3,L},27={R1,R3/2-R4,2*R3+R5,T1,T3,L},28={2*R1-2*R3+R4,2*R3+R5,T1,T3,L},29={R1,R3/2-R4,2*R3+R5,T1,L},30={R1,2*R3-R4,2*R3+R5,T1,L},31={R1,2*R3-R4,2*R3+R5,T3,L},32={R1,R4+R5,T1,T3,L},33={R1,R3/2-R4,R3/2+R5,T1,L}]):
*)
unitarity_constraints:=table([1={T2,-R2,-R3*(2*R3+R5)*(R3+2*R5)},2={T2,-R2,-R3*(2*R3+R5)*(R3+2*R5)},3={-R1*(R1+R5)*(2*R1+R5)},4={-R1*(R1+R5)*(2*R1+R5)},5={-R3*(2*R3+R5)*(R3+2*R5)},6={-R3*(2*R3+R5)*(R3+2*R5)},7={-R3*(2*R3+R5)*(R3+2*R5)},8={R1*(R1-2*R3-R5)*(2*R3+R5)},9={-R3*(2*R3+R5)*(R3+2*R5)},10={-R3*(2*R3+R5)*(R3+2*R5)},11={-R3*(2*R3+R5)*(R3+2*R5)},12={T2,-R2},13={-R2,-T1},14={T2,-R2},15={-R2,-T1},16={T2,-R2},17={T2,-R2},18={T2,-R2},19={T2,-R2},20={T2,-R2},21={T2,-R2},22={T2,-R2},23={T2,-R2},24={T2,-R2},25={-R2,-T1},26={T2,-R2},27={T2,-R2},28={T2,-R2},29={T2,-R2},30={T2,-R2},31={-R2,-T1},32={T2,-R2},33={T2,-R2}]);
#*)

cosmic_coefficients:=table([L="l",T1="t_1",T2="t_2",T3="t_3",R1="r_1",R2="r_2",R3="r_3",R4="r_4",R5="r_5",R6="r_6"]);

for jj from 1 to 33 do
  constraint:=eval(parameter_constraints[jj]);
  constraint:=map(primpart,constraint);
  consr:=convert(constraint,string);
  for ii in L,T1,T2,T3,R1,R2,R3,R4,R5,R6 do
    str:=eval(cosmic_coefficients[ii]);
    consr:=Substitute(consr,convert(ii,string),str);
    consr:=Substitute(consr,convert(ii,string),str);
  end do;
  star:=proc(s::string)::boolean;
    return evalb(s="*");
  end proc;
  consr:=Remove(star,consr);
  consr:=SubstituteAll(consr,",","=");
  consr:=cat("${",consr,"=0}$");
  constraint_||jj:=consr;

  constraint:=eval(unitarity_constraints[jj]);
  inequalitise:=proc(expr::algebraic)::`<`;
    if op(1,expr)=-1 then
       return -expr<0;	
    else
      return expr>0;
    end if;
  end proc;
  constraint:=map(inequalitise,constraint);
  consr:=convert(constraint,string);
  for ii in L,T1,T2,T3,R1,R2,R3,R4,R5,R6 do
    str:=eval(cosmic_coefficients[ii]);
    consr:=Substitute(consr,convert(ii,string),str);
    consr:=Substitute(consr,convert(ii,string),str);
    consr:=Substitute(consr,convert(ii,string),str);
  end do;
  star:=proc(s::string)::boolean;
    return evalb(s="*");
  end proc;
  consr:=Remove(star,consr);
  consr:=cat("${",consr,"}$");
  unitarity_constraint_||jj:=consr;
  print(%);
end do;
#fin();
particles:=proc(critical_case::integer)::NULL;
  local tmp,tmp2,gauge,gauges,m,p;
  tmp:=particle_content[critical_case];
  for JP in {J0Pm,J0Pp,J1Pm,J1Pp,J2Pm,J2Pp} do
    tmp2:=tmp[JP];
    if assigned(tmp2) then
      texfilename:=cat("/home/williamb/Documents/physics/papers/paper-2/paper/aps/particle_content_icons/critical_case_",convert(critical_case,string),"_",convert(JP,string),".tex");
      gauge:=-1;
      gauges:=numelems(convert(eval(tmp2),list));
      if gauges=1 then
	m:=3.5;
	p:=3.5;
      elif gauges=2 then
	m:=2.5;
	p:=4.5;
      elif gauges=3 then
	m:=1.5;
	p:=5.5;
      end if;
      m:=convert(m,string);
      p:=convert(p,string);
      writeto(texfilename);
      printf("\\documentclass[preview]{standalone}\n");
      printf("\\usepackage{tikz}\n");
      printf("\\usetikzlibrary{calc}\n");
      printf("\\usepackage{xcolor}\n");
      printf("\\definecolor{tordion}{RGB}{0,0,255}\n");
      printf("\\definecolor{agraviton}{RGB}{255,0,0}\n");
      printf("\\definecolor{sgraviton}{RGB}{0,255,0}\n");
      printf("\\begin{document}\n");
      printf("\\begin{tikzpicture}[node distance = 1.4cm, auto]\n");
      printf(cat("\\clip (225:",m,") rectangle (45:",p,");\n"));
      for excitation in tmp2 do
	gauge:=gauge+1;
	print_excitation(excitation,gauge);
      end do;
      printf("\\end{tikzpicture}\n");
      printf("\\end{document}\n");
      writeto(terminal);
    end if;
  end do;
  return NULL;
end proc;

sample_particles:=proc(critical_case::integer)::NULL;
  local tmp,tmp2,gauge,gauges,m,p;
  tmp:=particle_content[critical_case];
  for JP in {J0Pm,J0Pp,J1Pm,J1Pp,J2Pm,J2Pp} do
    tmp2:=tmp[JP];
    if assigned(tmp2) then
      texfilename:=cat("/home/williamb/Documents/physics/papers/paper-2/paper/aps/particle_content_icons/sample_critical_case_",convert(critical_case,string),"_",convert(JP,string),".tex");
      gauge:=-1;
      gauges:=numelems(convert(eval(tmp2),list));
      if gauges=1 then
	m:=2;
	p:=2;
      elif gauges=2 then
	m:=1;
	p:=3;
      elif gauges=3 then
	m:=0;
	p:=4;
      end if;
      m:=convert(m,string);
      p:=convert(p,string);
      writeto(texfilename);
      printf("\\documentclass[preview]{standalone}\n");
      printf("\\usepackage{tikz}\n");
      printf("\\usetikzlibrary{calc}\n");
      printf("\\usepackage{xcolor}\n");
      printf("\\definecolor{tordion}{RGB}{0,0,255}\n");
      printf("\\definecolor{agraviton}{RGB}{255,0,0}\n");
      printf("\\definecolor{sgraviton}{RGB}{0,255,0}\n");
      printf("\\begin{document}\n");
      printf("\\begin{tikzpicture}[node distance = 1.4cm, auto]\n");
      printf(cat("\\clip (225:",m,") rectangle (45:",p,");\n"));
      for excitation in tmp2 do
	gauge:=gauge+1;
	print_excitation(excitation,gauge);
      end do;
      printf("\\end{tikzpicture}\n");
      printf("\\end{document}\n");
      writeto(terminal);
    end if;
  end do;
  return NULL;
end proc;

print_excitation:=proc(excitation,gauge::integer)::NULL;
    local colour,x,y,r,coordinate,excitation_string;
    x:=convert(1.5*gauge,string);
    y:=convert(1.5*gauge,string);
    r:=convert(2*gauge,string);
    coordinate:=cat("gauge",convert(gauge,string));
    #printf(cat("\\coordinate (",coordinate,") at (",x,",",y,");\n"));
    printf(cat("\\coordinate (",coordinate,") at (45:",r,");\n"));
    if type(excitation,list) then
      excitation_string:=convert(excitation[1],string);
      colour:=colour_selector(excitation_string);
      if Search("v",excitation_string)<>0 then
	printf(cat("\\fill[",colour,"] (",coordinate,") -- (45:1) arc (45:225:1) -- cycle;\n"));
      elif Search("l",excitation_string)<>0 then
	printf(cat("\\fill[",colour,"] ($(",coordinate,")+(45:1)$) arc (45:225:1) -- ($(",coordinate,")+(225:0.5)$) -- ($(",coordinate,")+(225:0.5)$) arc (225:45:0.5) -- cycle;\n"));
      end if;
      excitation_string:=convert(excitation[2],string);
      colour:=colour_selector(excitation_string);
      if Search("v",excitation_string)<>0 then
	printf(cat("\\fill[",colour,"] (",coordinate,") -- (-135:1) arc (-135:45:1) -- cycle;\n"));
      elif Search("l",excitation_string)<>0 then
	printf(cat("\\fill[",colour,"] ($(",coordinate,")+(-135:1)$) arc (-135:45:1) -- ($(",coordinate,")+(45:0.5)$) -- ($(",coordinate,")+(45:0.5)$) arc (45:-135:0.5) -- cycle;\n"));
      end if;
    else
      excitation_string:=convert(excitation,string);
      colour:=colour_selector(excitation_string);
      if Search("v",excitation_string)<>0 then
	printf(cat("\\fill[",colour,"] (",coordinate,") circle (1);\n"));
      elif Search("l",excitation_string)<>0 then
	printf(cat("\\fill[",colour,",even odd rule] (",coordinate,") circle (1) (",coordinate,") circle (0.5);\n"));
      end if;
    end if;
  return NULL;
end proc;

colour_selector:=proc(excitation_string::string)::string;
  local colour;
  if Search("A",excitation_string)<>0 then
    colour:="tordion";
  elif Search("a",excitation_string)<>0 then
    colour:="agraviton";
  elif Search("s",excitation_string)<>0 then
    colour:="sgraviton";
  end if;
  return colour;
end proc;

for ii from 1 to 35 do
    if assigned(particle_content[ii]) then
     particles(ii);
     sample_particles(ii);
    end if;
end do;

#	now we would like to craft the flesh of the table

writeto("/home/williamb/Documents/physics/papers/paper-2/paper/aps/particle_content_table.tex");
for ii from 1 to 33 do
  if numelems(eval(case_translation[ii]))=1 then
    printf(cat(convert(case_translation[ii][1],string)," &",constraint_||ii));
  else
    prop1:=cat("$^{*",convert(case_translation[ii][2],string),"}$"):
    printf(cat(prop1,convert(case_translation[ii][1],string)," &",constraint_||ii));
  end if:
  printf(cat("&",unitarity_constraint_||ii));
  printf(cat("& \\includegraphics[width=0.4cm]{particle_content_icons/critical_case_",convert(ii,string),"_J0Pm.pdf}"));
  printf(cat("& \\includegraphics[width=0.4cm]{particle_content_icons/critical_case_",convert(ii,string),"_J0Pp.pdf}"));
  printf(cat("& \\includegraphics[width=0.4cm]{particle_content_icons/critical_case_",convert(ii,string),"_J1Pm.pdf}"));
  printf(cat("& \\includegraphics[width=0.4cm]{particle_content_icons/critical_case_",convert(ii,string),"_J1Pp.pdf}"));
  printf(cat("& \\includegraphics[width=0.4cm]{particle_content_icons/critical_case_",convert(ii,string),"_J2Pm.pdf}"));
  printf(cat("& \\includegraphics[width=0.4cm]{particle_content_icons/critical_case_",convert(ii,string),"_J2Pp.pdf}"));
  if ii=1 then
     printf(cat("& \\multirow{2}{*}{\\includegraphics[width=0.4cm]{massive.pdf}\\includegraphics[width=0.4cm]{massless.pdf}\\includegraphics[width=0.4cm]{massless.pdf}}"));
  elif ii=3 then
     printf(cat("& \\multirow{9}{*}{\\includegraphics[width=0.4cm]{massless.pdf}\\includegraphics[width=0.4cm]{massless.pdf}}"));
  elif ii=12 then
     printf(cat("& \\multirow{22}{*}{\\includegraphics[width=0.4cm]{massive.pdf}}"));
  else
    printf("&");
  end if;
  printf("\\\\ \n");
  if ii in {2,11} then
    printf("\\hline \n");
  end if;
end do;
writeto(terminal);

fin();












saved_particle_content:=table([1=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]]),2=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]]),3=table([J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[Al2],J2Pm=[Al2]]),4=table([J1Pm=[Al2],J1Pp=[Al2],J2Pm=[Al2]]),5=table([J1Pm=[AL2],J1Pp=[Al2],J2Pp=[Al2]]),6=table([J0Pm=[Al2],J1Pm=[Al2],J1Pp=[Al2],J2Pp=[Al2]]),7=table([J0Pm=[Al0],J1Pm=[Al2],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]]),8=table([J0Pp=[Al2],J1Pm=[Al2],J1Pp=[Al2],J2Pm=[Al2]]),9=table([J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[Al2],J2Pp=[Al2]]),10=table([J0Pm=[Al0],J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]]),11=table([J0Pm=[Al2],J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[Al2],J2Pp=[Al2]]),12=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al0,Al0},{Al0,sl2},{Al0,al2}],J1Pp=[{Al0,Al0},{Al0,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),13=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al0,Al0},{Al0,sl2},{Al0,al2}],J1Pp=[{Ali,Al0},{Ali,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),14=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Ali,Al0},{Ali,sl2},{Ali,al2}],J1Pp=[{Al0,Al0},{Al0,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),15=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al0,Al0},{Al0,sl2},{Al0,al2}],J1Pp=[{Al0,Al0},{Al0,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),16=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[{Al2,Al0},{Al2,al2}]]),17=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[Al0,sl2,al2],J1Pp=[Al0,al2]]),18=table([J0Pm=[Av2],J1Pp=[Al0,al2]]),19=table([J0Pm=[Av2],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]]),20=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[{Al2,Al0},{Al2,al2}]]),21=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[Al0,al2],J2Pm=[Al2]]),22=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[Al0,al2],J2Pm=[Al2]]),23=table([J0Pm=[Av2],J0Pp=[Al2],J1Pp=[Al0,al2]]),24=table([J0Pm=[Av2],J1Pm=[Al2,sl2,al2],J1Pp=[{Al0,al2},{Al0,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),25=table([J0Pm=[Av2],J1Pm=[Al2,sl2,al2],J1Pp=[{Ali,al2},{Ali,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),26=table([J0pm=[Av2],J0Pp=[Al2],J1Pm=[Al2],J1Pp=[{Al2,Al0},{Al2,al2}]]),27=table([J0Pm=[Av2],J1Pm=[Al2],J1Pp=[Al0,al2],J2Pp=[Al2]]),28=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al2],J1Pp=[Al0,al2],J2Pm=[Al2]]),29=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[{Al2,Al0},{Al2,sl2},{Al2,al2}],J1Pp=[Al0,al2],J2Pp=[Al2]]),30=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al0,sl2,al2],J1Pp=[{Al0,Al0},{Al0,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),31=table([J0Pm=[Av2],J0Pp=[Al2],J1Pm=[Al0,sl2,al2],J1Pp=[{Ali,Al0},{Ali,al2}],J2Pm=[Al0],J2Pp=[Al0,sl2]]),32=table([J0Pm=[Av2],J0Pp=[Al2],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]]),33=table([J0Pm=[Av2],J0Pp=[Al0,sl2],J1Pm=[Al0,sl2,al2],J1Pp=[{Al2,Al0},{Al2,al2}],J2Pp=[Al2]])]);
