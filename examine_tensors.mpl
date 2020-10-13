
restart:
with(Physics):
read `sta.mpl`:
with(sta):

#	alias the dynamical fields
alias(a=S(t),b=conformal_isometry(t),V=dilation_gauge(t),X=normal_torsion(t),Y=dual_torsion(t),Phi=compensator(t)):
#	define isotropic torsion
Z:=ds((X+diff(a*b,t)/b)/(a*b)+(Y/(a*b))*ps):
#	define gravity frames
gu0:=ds(eu0/b):
gu1:=ds(sqrt(1-k*r^2)/(a*b)*eu1):
gu2:=ds(1/(a*b)*eu2):
gu3:=ds(1/(a*b)*eu3):
dirac[make_gd]():
#	define the dilation gauge field
calV_:=ds(V/(a*b)*gd0):
V:=0:

read `dirac/dirac_maple_general_h_functions_v6.mpl`:
read `ewgt.mpl`:
with(ewgt):
read `tools.mpl`:
with(tools):

ewgt[form_eWGT_quantities]():
ewgt[form_quadratic_Riemann]():
ewgt[form_quadratic_torsion]():
ewgt[form_miscellaneous]():


tools[print_for_paper]():
