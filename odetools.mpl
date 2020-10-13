#	this module get things out of the way for the numerical

Erato:=module() 
description "some stuff for numerical";
option package;
export convert_to_conformal, convert_to_simple_conformal, convert_to_dimensionless, render_inert, py_syntax;

  convert_to_conformal:=proc(expr)
    description "converts any expression to conformal time";
    local tmp2;
    if type(expr,equation) then
      tmp2:=dchange(t=tfun(x),expr,{x},known={tfun});
      tmp2:=evala(simplify(subs(diff(tfun(x),x)=a(x)/(s*d0a),tmp2)));
      return simplify(tmp2);
    else
      tmp2:=dchange(t=tfun(x),expr,{x},known={tfun});
      tmp2:=evala(simplify(subs(diff(tfun(x),x)=a(x)/(s*d0a),tmp2)));
      return simplify(numer(tmp2));
    end if;
  end proc:

  convert_to_simple_conformal:=proc(expr)
    description "converts any expression to trivial conformal time";
    local tmp2;
    if type(expr,equation) then
      tmp2:=dchange(t=tfun(x),expr,{x},known={tfun});
      tmp2:=evala(simplify(subs(diff(tfun(x),x)=a(x)/(S*M_p),tmp2)));
      return simplify(tmp2);
    else
      tmp2:=dchange(t=tfun(x),expr,{x},known={tfun});
      tmp2:=evala(simplify(subs(diff(tfun(x),x)=a(x)/(S*M_p),tmp2)));
      return simplify(numer(tmp2));
    end if;
  end proc:

  render_inert:=proc(expr)
    description "takes an expression and converts to contemporary inert form";
    local tmp1;
    tmp1:=simplify(subs(diff(diff(a(x),x),x)=d2a,diff(diff(Q(x),x),x)=d2Q,expr));
    tmp1:=simplify(subs(diff(a(x),x)=d1a,diff(Q(x),x)=d1Q,tmp1));
    tmp1:=simplify(subs(R(x)=R0,H(x)=H0,q(x)=q0,O_r(x)=O_r0,O_d(x)=O_d0,O_L(x)=O_L0,a(x)=d0a,Q(x)=d0Q,tmp1));
    return tmp1;
  end proc;

  convert_to_dimensionless:=proc(expr)
    description "takes an expression and converts K=l^2 (l is length scale), then s=l*H0 (s is dimensionless scaling of conformal time which picks out conformal Hubble times), and finally l=R0=1";
    local tmp1;
    if type(expr,equation) then
      tmp1:=simplify(subs(K=l^2,(expr))):		#replace K by the length unit
      tmp1:=simplify(subs(s=l*H0,tmp1)):	#replace s by the unit that switches to conf hubble times
      return simplify(subs(R0=1,l=1,tmp1)):		#set length unit to 1
    else
      tmp1:=simplify(subs(K=l^2,(expr))):		#replace K by the length unit
      tmp1:=simplify(subs(s=l*H0,tmp1/(H0^2))):	#replace s by the unit that switches to conf hubble times
      return simplify(subs(R0=1,l=1,tmp1)):		#set length unit to 1
    end if:
  end proc:

  py_syntax:=proc(object)::string;
    description "takes a maple object and converts the syntax to python";
    local temp;
    writeto("tmp.mpl");
    lprint(eval(object));
    print(`;`);
    writeto(terminal);
    read(`tmp.mpl`);
    temp:=SubstituteAll(convert(%,string),"^","**");
    temp:=SubstituteAll(temp,"[1]","[0]");
    temp:=SubstituteAll(temp,"[2]","[1]");
    temp:=SubstituteAll(temp,"[3]","[2]");
    temp:=SubstituteAll(temp,"[4]","[3]");
    temp:=SubstituteAll(temp,"[5]","[4]");
    temp:=SubstituteAll(temp,"[6]","[5]");
    return temp;
  end proc;

end module;
