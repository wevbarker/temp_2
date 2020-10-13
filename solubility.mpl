f(x):=3*x^2+10*x+3:
prob:=(A1+A2*f(x)*diff(f(x),x))^2-A3*f(x)*diff(f(x),x)^2=0:
collect(expand(prob),x);
simplify(%);
solve(%,x);

