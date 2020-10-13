read(`theory_tools.mpl`):

for ii in Y,M,W,Z do
  for jj in Y,M,W,Z do
    if evalb(ii<>jj) then
      convert_parameters_print_equations(ii,jj);
    end if:
  end do:
end do:

dit({},"here are constraints we have found imposed by various authors");

testcase_literature(Minkevich);
testcase_literature(Zhang);
testcase_literature(SNY1);
testcase_literature(SNY2);
testcase_literature(pure_weyl_GA);
testcase_literature(pure_weyl_GA_1);
testcase_literature(pure_weyl_GA_2);
testcase_literature(pure_weyl_GA_3);
convert_parameters_print_equations(cW1,M);

translate_cosmological_parameters(cW1,cZ);

fin();

(*
dit({},"D");
testcase_critical_case(16);
dit({},"H");
testcase_critical_case(14);
dit({},"J");
testcase_critical_case(11);
dit({},"O");
testcase_critical_case(10);
dit({},"E");
testcase_critical_case(1);
dit({},"E extras");
testcase_critical_case(27);
testcase_critical_case(30);
testcase_critical_case(35);
dit({},"of interest...");
testcase_critical_case(15);
testcase_critical_case(12);
*)

dit({},"here are cosmological coordinate translations");

translate_cosmological_parameters(cW1,cV);
translate_cosmological_parameters(cW1,cZ);
convert_parameters_print_equations(cW1,G);


convert_parameters_print_equations(M,Y);
convert_parameters_print_equations(Y,M);
convert_parameters_print_equations(M,NY);
convert_parameters_print_equations(NY,M);


fin();
