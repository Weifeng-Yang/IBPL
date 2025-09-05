%% Input.
% var         : Initial matrix
% ngmar       : Decomposed tensor
% aa          : Maximum number of non-zero elements for each decomposition matrix
% The remaining parameters are explained the same as the 'main_Run_me' function

%% Output.
% vars:       : Decomposition matrix resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         : Runtime array during iteration
% btss and atss: An array of extrapolation parameters produced by each algorithm during iteration,


%% An Inertial Block Proximal Linearized Method with Adaptive Momentum for Nonconvex and Nonsmooth Optimization

function [data,varss]=ALGOchoose(var,ngmar,aa,maxiteropt,at,bt,flag,stopindex,r,t,btmax,rho)

if(flag==1) 
[vars,loss,tr,bts,ats]=IPALM(var,ngmar,aa,maxiteropt,stopindex);
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;
atss=ats;

elseif(flag==2)
[vars,loss,tr,bts]=BPL(var,ngmar,aa,maxiteropt,stopindex,r);
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;


elseif(flag==3) 
[vars,loss,tr,bts,ats]=IBPG(var,ngmar,aa,maxiteropt,stopindex);
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;
atss=ats;

elseif(flag==4) 
[vars,loss,tr,bts]=TITAN(var,ngmar,aa,maxiteropt,stopindex);
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;


elseif(flag==5)
[vars,loss,tr,bts]=ABPL(var,ngmar,aa,maxiteropt,bt,1,stopindex,t,btmax);  
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;

elseif(flag==6)
[vars,loss,tr,bts]=PGels(var,ngmar,aa,maxiteropt,stopindex,r);  
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;

elseif(flag==7) 
[vars,loss,tr,bts]=APGL(var,ngmar,aa,maxiteropt,stopindex,btmax,rho);  
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;

elseif(flag==8)
[vars,loss,tr,ats,bts]=IBPL(var,ngmar,aa,maxiteropt,at,bt,2,stopindex,r,t,btmax,rho);  %%IBPL
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;
atss=ats;


elseif(flag==9)
[vars,loss,tr,ats,bts]=IBPL(var,ngmar,aa,maxiteropt,at,bt,1,stopindex,r,t,btmax,rho);  %%IBPL+
varss=vars;
lossdata=loss;
trdata=tr;
btss=bts;
atss=ats;

end


data{1}=lossdata;
data{2}=trdata;
data{3}=btss;
end