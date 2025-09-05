%% The parameters update function of IBPL
function var=IBPLupdate(vv,varz,var,varze,num,ngmar,aa,r)
for j=1:num
V=gradibpl(vv,varz,var,ngmar,j,num,r);
[var{j},~]=PROX(varze{j},V,aa(j));
end
end