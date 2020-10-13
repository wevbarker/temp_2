#	this module get things out of the way for the numerical

Erato:=module() 
description "some stuff for numerical";
option package;
export convert_to_conformal, render_inert, py_syntax, wl_syntax, wl_syntax_list;

  convert_to_conformal:=proc(expr)
    description "converts any expression to conformal time";
    local tmp2;
    if type(expr,equation) then
      tmp2:=dchange(t=tfun(x),expr,{x},known={tfun});
      tmp2:=evala(simplify(subs(diff(tfun(x),x)=a(x),tmp2)));
      return simplify(tmp2);
    else
      tmp2:=dchange(t=tfun(x),expr,{x},known={tfun});
      tmp2:=evala(simplify(subs(diff(tfun(x),x)=a(x),tmp2)));
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

  py_syntax:=proc(object)::string;
    description "takes a maple object and converts the syntax to python";
    writeto("tmp.mpl");
    lprint(eval(object));
    print(`;`);
    writeto(terminal);
    read(`tmp.mpl`);
    temp:=SubstituteAll(convert(%,string),"^","**");
    temp:=SubstituteAll(temp,"cosh(","np.cosh(");
    temp:=SubstituteAll(temp,"arccosh(","np.arccosh(");
    temp:=SubstituteAll(temp,"sinh(","np.sinh(");
    temp:=SubstituteAll(temp,"arcsinh(","np.arcsinh(");
    temp:=SubstituteAll(temp,"[1]","[0]");
    temp:=SubstituteAll(temp,"[2]","[1]");
    temp:=SubstituteAll(temp,"[3]","[2]");
    temp:=SubstituteAll(temp,"[4]","[3]");
    return temp;
  end proc;

  wl_syntax:=proc(object)::string;
    description "takes a maple object and converts the syntax to wolfram language";
    icosh:=proc(expr)
      return Cosh[expr]:
    end proc:
    isinh:=proc(expr)
      return Sinh[expr]:
    end proc:
    iarccosh:=proc(expr)
      return ArcCosh[expr]:
    end proc:
    iarcsinh:=proc(expr)
      return ArcSinh[expr]:
    end proc:
    writeto("tmp.mpl");
    lprint(eval(object));
    print(`;`);
    writeto(terminal);
    read(`tmp.mpl`);
    temp:=simplify(subs(sinh=isinh,cosh=icosh,arcsinh=iarcsinh,arccosh=iarccosh,%));
    temp:=map(simplify,temp);
    temp:=map(eval,temp);
    temp:=convert(temp,string);
    temp:=SubstituteAll(temp,"=","==");
    temp:=SubstituteAll(temp,"Y[1]","Y1[u]");
    temp:=SubstituteAll(temp,"Y[2]","Y2[u]");
    temp:=SubstituteAll(temp,"Y[3]","Y3[u]");
    temp:=SubstituteAll(temp,"Y[4]","Y4[u]");
    temp:=SubstituteAll(temp,"Y[5]","Y5[u]");
    temp:=SubstituteAll(temp,"YP[1]","Y1'[u]");
    temp:=SubstituteAll(temp,"YP[2]","Y2'[u]");
    temp:=SubstituteAll(temp,"YP[3]","Y3'[u]");
    temp:=SubstituteAll(temp,"YP[4]","Y4'[u]");
    temp:=SubstituteAll(temp,"YP[5]","Y5'[u]");
    return temp;
  end proc;

  wl_syntax_list:=proc(object)::string;
    description "takes a maple object and converts the syntax to wolfram language";
    writeto("tmp.mpl");
    lprint(eval(object));
    print(`;`);
    writeto(terminal);
    read(`tmp.mpl`);
    temp:=convert(%,string);
    temp:=SubstituteAll(temp,"[1]","[0]");
    temp:=SubstituteAll(temp,"[2]","[1]");
    temp:=SubstituteAll(temp,"[3]","[2]");
    temp:=SubstituteAll(temp,"[4]","[3]");
    temp:=SubstituteAll(temp,"[","{");
    temp:=SubstituteAll(temp,"]","}");
    return temp;
  end proc;

end module;
