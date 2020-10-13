#	The general-purpose main function

restart;

with(Physics):
with(PDEtools):
with(CLIo);
read `dirac/dirac_maple_test_r6.mpl`:
read `dirac/dirac_maple_general_h_functions_v6.mpl`:
read `dirac/spherical_gravity_for_r10.mpl`:
read `sta.mpl`:
with(sta):
read `ewgt.mpl`:
with(ewgt):
read `tools.mpl`:
with(tools):

case:=PGT_case_3;
theory[gauge_theory]:=PGT;
theory[Ricci]:=A;
theory[quadratic_Riemann]:=S1,S2,S3;
theory[quadratic_torsion]:=U1,U2;
theory[dust]:=Rho_d;
theory[radiation]:=Rho_r;
theory[cosmological_constant]:=Ls;
theory[curvature]:=k;

tools[refine_equations_of_motion]();
#tools[walkthrough](analytic);
