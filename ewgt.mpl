#	the procedures required for eWGT calculations

#	the ebui and ebdi are Lorentz basis vectors

ebd0:= that; ebd1:= rhat; ebd2:= thetahat; ebd3:=phihat;
ebu0:= that; ebu1:= rhat; ebu2:= thetahat; ebu3:=phihat;

ewgt[form_metrics]:=proc()
local i,j;
  for i from 0 to 3 do
    for j from 0 to 3 do
      ed_int_ed_||i||j:=ds(ed||i&@ed||j)[1];
      eu_int_eu_||i||j:=ds(eu||i&@eu||j)[1];
      gd_int_gd_||i||j:=ds(gd||i&@gd||j)[1];
      gu_int_gu_||i||j:=ds(gu||i&@gu||j)[1];
    od;
  od;
end:

#	some bulk procedures

ewgt[form_eWGT_quantities]:=proc()
  description "forms many quantities which will be of interest in eWGT";
  comment({},"forming eWGT quantities...");
  ewgt[form_metrics]():
  ewgt[form_calT_ed_]():
  ewgt[form_calT_gd_]():
  ewgt[form_calT_]():
  ewgt[form_T_]():
  ewgt[form_eui_ext_calT_edi_]():
  ewgt[form_Omi]():
  ewgt[form_omi]():
  ewgt[form_V_]():
  ewgt[form_primed_Omi]():
  ewgt[form_R_dagger_ed_ed_]():
  ewgt[form_calR_dagger_ed_ed_]():
  ewgt[form_eui_int_calR_dagger_edi_ed_]():
  ewgt[form_eui_ext_calR_dagger_edi_ed_]():
  ewgt[form_euj_int_eui_int_calR_dagger_edi_edj_]():
  ewgt[form_eui_int_calR_dagger_bar_edi_ed_]():
  ewgt[form_calR_dagger_bar_ed_ed_]():
  ewgt[form_calT_dagger_ed_]():
  ewgt[form_eui_ext_calT_dagger_edi_]():
  ewgt[form_calT_dagger_]():
  ewgt[form_calD_dagger_phi_]():
  ewgt[form_testify_edi_]():
  comment({},"...done");
end proc;
 
ewgt[form_quadratic_Riemann]:=proc()
  description "forms the contributions to the eWGT lagrangian quadratic in the Riemann tensor";
  comment({},"forming quadratic Riemann invariants...");
  ewgt[form_einstein_hilbert_]():
  ewgt[form_alpha_1_]():
  ewgt[form_alpha_2_]():
  ewgt[form_alpha_3_]():
  ewgt[form_alpha_4_]():
  ewgt[form_alpha_5_]():
  ewgt[form_alpha_6_]():
  comment({},"...done");
end proc;

ewgt[form_quadratic_torsion]:=proc()
  description "forms the contributions to the eWGT lagrangian quadratic in the torsion tensor";
  comment({},"forming quadratic torsion invariants...");
  ewgt[form_beta_1_]():
  ewgt[form_beta_2_]():
  ewgt[form_beta_3_]():
  ewgt[form_beta_dagger_1_]():
  ewgt[form_beta_dagger_2_]():
  ewgt[form_beta_dagger_3_]():
  comment({},"...done");
end proc;

ewgt[form_miscellaneous]:=proc()
  description "forms the contributions to the eWGT lagrangian from the compensator and matter";
  comment({},"forming miscellaneous invariants...");
  ewgt[form_lambda_]():
  ewgt[form_nu_]():
  ewgt[form_radiation_lagrangian]():
  ewgt[form_dust_lagrangian]():
  comment({},"...done");
end proc;

#	basic quantities in eWGT and PGT

#       torsion quantities

ewgt[form_calT_ed_]:=proc()
        local ii, tmp;
        for ii from 0 to 3 do
                tmp:=dbpart(ed||ii&@ed0);
                calT_ed_||ii:=dbpart(Z&@tmp);
                tmp:=dbpart(gammau||ii&@ed0);
                calT_gammau_||ii:=dbpart(Z&@tmp);
        od;
end:

ewgt[form_calT_gd_]:=proc()
  local ii, jj, tmp;
  for ii from 0 to 3 do
    tmp:=ds(0);
    for jj from 0 to 3 do
      tmp:=tmp+ds(ds(gd||ii&@eu||jj)[1]*calT_ed_||jj);
    od;
    calT_gd_||ii:=tmp;
  od;
end:

ewgt[form_calT_]:=proc()
  local tmp, jj;
  global calT_;
  tmp:=ds(0);
  for jj from 0 to 3 do
    tmp:=tmp+dvpart(eu||jj&@calT_ed_||jj);
  od;
  calT_:=ds(tmp);
end:

ewgt[form_eui_ext_calT_edi_]:=proc()
  local tmp,ii;
  global eui_ext_calT_edi_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    tmp:=tmp+dtpart(eu||ii&@calT_ed_||ii);
  od;
  eui_ext_calT_edi_:=ds(tmp);
end proc;

ewgt[form_testify_edi_]:=proc()
  description "forms a test quantity, please ignore this!":
  global Acoef,Bcoef;
  local ii;
  print("here are the four bivector components of testify:");
  Acoef:=ds(A):
  Bcoef:=ds(A):
  for ii from 0 to 3 do
    testify_ed||ii||_:=dbpart(A&@(calT_ed_||ii+calT_&@ed||ii)+B&@(calT_ed_||ii-ed||ii&@eui_ext_calT_edi_)):
    print(ds(testify_ed||ii||_));
    #print("and calT_&@ed_||ii is");
    #print(ds(calT_&@(ed_||ii)));
  end do:
  return NULL:
end proc:








ewgt[form_T_]:=proc()
  local tmp, ii;
  global T_;
  tmp:=0;
  for ii from 0 to 3 do
    tmp:= ds(tmp + eu||ii*ds(gd||ii&@calT_)[1]);
  od;
  T_:=tmp;
end:

ewgt[form_calT_dagger_ed_]:=proc()
  local ii;
  for ii from 0 to 3 do
    calT_dagger_ed_||ii:=ds(calT_ed_||ii+(1/3)*dbpart(calT_&@ed||ii));
    #print(%);
  od;
end:

ewgt[form_eui_ext_calT_dagger_edi_]:=proc()
  local tmp,ii;
  global eui_ext_calT_dagger_edi_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    tmp:=tmp+dtpart(eu||ii&@calT_dagger_ed_||ii);
    #print(%);
  od;
  eui_ext_calT_dagger_edi_:=ds(tmp);
  #print(ds(tmp));
end proc;

ewgt[form_calT_dagger_]:=proc()
  local tmp, jj;
  global calT_dagger_;
  tmp:=ds(0);
  for jj from 0 to 3 do
    tmp:=tmp+dvpart(eu||jj&@calT_dagger_ed_||jj);
  od;
  calT_dagger_:=ds(tmp);
  #print(%);
end:

#	some miscellaneous quantities

ewgt[form_V_]:=proc()
  local tmp, ii;
  global V_;
  tmp:=0;
  for ii from 0 to 3 do
    tmp:= ds(tmp + eu||ii*ds(gd||ii&@calV_)[1]);
  od;
  V_:=tmp;
end:

ewgt[form_calH_dagger_]:=proc()
  local tmp, ii;
  global calH_;
  tmp:=0;
  for ii from 0 to 3 do
    tmp:=tmp-eu||ii&@map(diff,(V_+(1/3)*U_),x||ii);
  od;
  calH_:=hob(dbpart(tmp));
end:

ewgt[form_calD_dagger_phi_]:=proc()
  local tmp, ii;
  global calD_dagger_phi_;
  tmp:=0;
  for ii from 0 to 3 do
    tmp:=ds(tmp+gu||ii&@diff(compensator(t),x||ii));
  od;
  calD_dagger_phi_:=ds(tmp+compensator(t)*(V_+(1/3)*calT_));
end:

#       procedure to form the rotational gauge fields from the displacement gauge fields and the torsion

ewgt[form_Omi] := proc()

        local i,j,tmp,dot,tmp1,tmp2;

        for j from 0 to 3 do
                tmp:=0;
                for i from 0 to 3 do
                        tmp:= tmp+gu||i&@map(diff,gu||j,x||i);
                od;
                Bv||j:=dbpart(tmp);
        od;
        for j from 0 to 3 do
                tmp:=ds(0);
                for i from 0 to 3 do
                        dot:=ds(gd||i&@gd||j);
                        tmp2:=dot[1];
                        tmp1:=dvpart(gd||j&@Bv||i);
                        tmp:= ds(tmp+1/2*gd||i&@tmp1+1/2*tmp2*Bv||i);
                od;
                hOm||j:=dbpart(tmp);
                Om||j:=ds(dbpart(hOm||j-calT_gd_||j+(1/2)*gd||j&@eui_ext_calT_edi_));
        od;
end:

#       procedure to form the \omega(e_\mu) fields

ewgt[form_omi]:=proc()
  local ii,jj,tmp,tmp1;
  for ii from 0 to 3 do
      tmp:=hub(ed||ii);
      om||ii:=ds(0);
      for jj from 0 to 3 do
        tmp1:=ds(tmp&@eu||jj)[1];
        om||ii:=ds(om||ii+tmp1*Om||jj);
      od;
  od;
end:

ewgt[form_primed_Omi]:=proc()
  local ii;
  for ii from 0 to 3 do
    primed_Om||ii:=dbpart(Om||ii+calV_&@gd||ii);
  od;
end:

#	forming the contractions and protractions of the Riemann

ewgt[form_R_dagger_ed_ed_] := proc()
  local i,j,tmp;
  for i from 0 to 3 do
    for j from 0 to 3 do
      tmp := map(diff,primed_Om||j,x||i)-map(diff,primed_Om||i,x||j)+1/2*(primed_Om||i&@primed_Om||j-primed_Om||j&@primed_Om||i);
      R_dagger_ed_ed_||i||j := ds(tmp);
    od;
  od;
end:

ewgt[form_calR_dagger_ed_ed_]:=proc()
  local i,j,k,l,tmp;
  for i from 0 to 3 do
    for j from 0 to 3 do
      tmp:=0;
      for k from 0 to 3 do
        for l from 0 to 3 do
          tmp:=ds(tmp+ds(ed||i&@gu||k)[1]*ds(ed||j&@gu||l)[1]*R_dagger_ed_ed_||k||l);
        od;
      od;
      calR_dagger_ed_ed_||i||j:=tmp;
    od;
  od;
end:

ewgt[form_eui_int_calR_dagger_edi_ed_]:=proc()
  local ii, jj, tmp;
  for jj from 0 to 3 do
    tmp:=ds(0);
    for ii from 0 to 3 do
      tmp:=tmp+dvpart(eu||ii&@calR_dagger_ed_ed_||ii||jj);
    od;
    eui_int_calR_dagger_edi_ed_||jj:=ds(tmp);
  od;
end:

ewgt[form_eui_ext_calR_dagger_edi_ed_]:=proc()
  local ii, jj, tmp;
  for jj from 0 to 3 do
    tmp:=ds(0);
    for ii from 0 to 3 do
      tmp:=tmp+dtpart(eu||ii&@calR_dagger_ed_ed_||ii||jj);
    od;
    eui_ext_calR_dagger_edi_ed_||jj:=ds(tmp);
  od;
end:

ewgt[form_euj_int_eui_int_calR_dagger_edi_edj_]:=proc()
  local jj, tmp;
  global euj_int_eui_int_calR_dagger_edi_edj_;
  tmp:=ds(0);
  for jj from 0 to 3 do
    tmp:=tmp+ds(eu||jj&@eui_int_calR_dagger_edi_ed_||jj)[1];
  od;
  euj_int_eui_int_calR_dagger_edi_edj_:=ds(tmp);
end:

ewgt[form_eui_int_calR_dagger_bar_edi_ed_]:=proc()
  local ii, jj, tmp;
  for ii from 0 to 3 do
    tmp:=ds(0);
    for jj from 0 to 3 do
      tmp:=tmp+ds(ed||ii&@eui_int_calR_dagger_edi_ed_||jj)[1]*eu||jj;
    od;
    eui_int_calR_dagger_bar_edi_ed_||ii:=ds(tmp);
  od;
end:

ewgt[form_calR_dagger_bar_ed_ed_]:=proc()
  local ii, jj, kk, ll, tmp;
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      tmp:=ds(0);
      for kk from 0 to 3 do
	for ll from 0 to 3 do
	  tmp:=tmp+(1/2)*ds(dbpart(ed||ii&@ed||jj)&@calR_dagger_ed_ed_||kk||ll)[1]*dbpart(eu||ll&@eu||kk);
	od;
      od;
      calR_dagger_bar_ed_ed_||ii||jj:=ds(tmp);
    od;
  od;
end:

#	contributions to the eWGT lagrangia
#	contributions of the quadratic Riemann

ewgt[form_einstein_hilbert_]:=proc()
  global einstein_hilbert_;
  einstein_hilbert_:=ds(euj_int_eui_int_calR_dagger_edi_edj_)[1];
end proc;

ewgt[form_alpha_1_]:=proc()
  global alpha_1_;
  alpha_1_:=ds(euj_int_eui_int_calR_dagger_edi_edj_&@euj_int_eui_int_calR_dagger_edi_edj_)[1];
end:

ewgt[form_alpha_2_]:=proc()
  local ii, jj, tmp;
  global alpha_2_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      tmp:=tmp+ds(eu_int_eu_||ii||jj*eui_int_calR_dagger_edi_ed_||ii&@eui_int_calR_dagger_edi_ed_||jj)[1];
    od;
  od;
  alpha_2_:=ds(tmp)[1];
end:

ewgt[form_alpha_3_]:=proc()
  local ii, jj, tmp, tmp2;
  global alpha_3_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      tmp:=tmp+ds(eu_int_eu_||ii||jj*eui_int_calR_dagger_bar_edi_ed_||ii&@eui_int_calR_dagger_edi_ed_||jj)[1];
    od;
  od;
  alpha_3_:=ds(tmp)[1];
end:

ewgt[form_alpha_4_]:=proc()
  local ii, jj, kk, ll, tmp;
  global alpha_4_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      for kk from 0 to 3 do
	for ll from 0 to 3 do
	  tmp:=tmp+ds(eu_int_eu_||ii||ll*eu_int_eu_||jj||kk*calR_dagger_ed_ed_||ii||jj&@calR_dagger_ed_ed_||kk||ll)[1];
	od;
      od;
    od;
  od;
  alpha_4_:=ds(tmp)[1];
end:

ewgt[form_alpha_5_]:=proc()
  local ii, jj, tmp;
  global alpha_5_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      tmp:=tmp+ds(eu_int_eu_||ii||jj*eui_ext_calR_dagger_edi_ed_||ii&@eui_ext_calR_dagger_edi_ed_||jj)[1];
    od;
  od;
  alpha_5_:=ds(tmp)[1];
end:

ewgt[form_alpha_6_]:=proc()
  local ii, jj, kk, ll, tmp;
  global alpha_6_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      for kk from 0 to 3 do
	for ll from 0 to 3 do
	  tmp:=tmp+ds(eu_int_eu_||ii||ll*eu_int_eu_||jj||kk*calR_dagger_bar_ed_ed_||ii||jj&@calR_dagger_ed_ed_||kk||ll)[1];
	od;
      od;
    od;
  od;
  alpha_6_:=ds(tmp)[1];
end:

#	terms quadratic in the torsion
#	eWGT torsion

ewgt[form_beta_dagger_1_]:=proc()
  local ii, jj, tmp;
  global beta_dagger_1_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      tmp:=tmp+ds(eu_int_eu_||ii||jj*calT_dagger_ed_||ii&@calT_dagger_ed_||jj)[1];
    od;
  od;
  beta_dagger_1_:=ds(tmp)[1];
end:

ewgt[form_beta_dagger_2_]:=proc()
  global beta_dagger_2_;
  beta_dagger_2_:=ds(eui_ext_calT_dagger_edi_&@eui_ext_calT_dagger_edi_)[1];
end proc;

ewgt[form_beta_dagger_3_]:=proc()
  global beta_dagger_3_;
  beta_dagger_3_:=ds(calT_dagger_&@calT_dagger_)[1];
end:

#	PGT torsion

ewgt[form_beta_1_]:=proc()
  local ii,jj,tmp;
  global beta_1_;
  tmp:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      tmp:=tmp+ds(eu_int_eu_||ii||jj*calT_ed_||ii&@calT_ed_||jj)[1];
    od;
  od;
  beta_1_:=ds(tmp)[1];
end proc;

ewgt[form_beta_2_]:=proc()
  global beta_2_;
  beta_2_:=ds(eui_ext_calT_edi_&@eui_ext_calT_edi_)[1];
end proc;

ewgt[form_beta_3_]:=proc()
  global beta_3_;
  beta_3_:=ds(calT_&@calT_)[1];
end proc;

#	the compensator kinetic term

ewgt[form_nu_]:=proc()
  global nu_;
  nu_:=ds(calD_dagger_phi_&@calD_dagger_phi_)[1];
end proc;

#	the compensator interaction term

ewgt[form_lambda_]:=proc()
  global lambda_;
  lambda_:=L;
end proc;

#	the radiation and the baryonic dust

ewgt[form_radiation_lagrangian]:=proc()
  global radiation_lagrangian;
  radiation_lagrangian:=-Rho_r/S(t)^4;
end proc;

ewgt[form_dust_lagrangian]:=proc()
  global dust_lagrangian;
  dust_lagrangian:=-Rho_d/S(t)^3;
end proc;

#	a procedure to find the Gauss-Bonnet divergence

machines[form_cal_gauss_bonnet_]:=proc()
  local ii,jj,kk,ll,tmp1,tmp2,tmp3;
  global cal_gauss_bonnet_;
  tmp1:=ds(0);
  tmp2:=ds(0);
  tmp3:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      for kk from 0 to 3 do
	for ll from 0 to 3 do
	tmp1:=ds(ps*ds(calR_dagger_ed_ed_||kk||ll&@calR_dagger_ed_ed_||ii||jj)[16]);
	tmp2:=ds(ps*ds(eu||ii&@eu||jj&@eu||kk&@eu||ll)[16]);
	tmp3:=tmp3+ds(tmp2&@tmp1)[1];
	od;
      od;
    od;
  od;
  cal_gauss_bonnet_:=ds(tmp3)[1];
end proc;

machines[form_gauss_bonnet_]:=proc()
  local ii,jj,kk,ll,tmp1,tmp2,tmp3,tmp4;
  global gauss_bonnet_;
  tmp1:=ds(0);
  tmp2:=ds(0);
  tmp3:=ds(0);
  tmp4:=ds(0);
  for ii from 0 to 3 do
    for jj from 0 to 3 do
      for kk from 0 to 3 do
	tmp1:=ds(ps*ds(R_dagger_ed_ed_||ii||kk&@Om||jj+(1/3)*Om||ii&@Om||jj&@Om||kk)[16]);
        tmp2:=ds(0);
	for ll from 0 to 3 do
	  tmp2:=tmp2+dtpart(eu||ll&@map(diff,tmp1,x||ll));
	od;
        tmp3:=dtpart(eu||ii&@eu||jj&@eu||kk);
        tmp4:=tmp4+ds(tmp3&@tmp2)[1];
      od;
    od;
  od;
  gauss_bonnet_:=ds(2*tmp4)[1];
end proc;

