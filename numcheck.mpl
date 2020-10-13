with(plots):
with(FileTools):
with(PDEtools):
with(DEtools):
with(StringTools):
with(SolveTools):
with(ODEsys);
read(`odetools.mpl`):
with(Erato);
with(CLIo):
read `tools.mpl`:
with(tools):
read `theory_tools.mpl`:

for ii from 1 to 1000 do
  jj:=2*ii/((4*ii-1)):
  if type(jj,integer) then
    print(ii,jj):
  end if:
end do:

fin();
