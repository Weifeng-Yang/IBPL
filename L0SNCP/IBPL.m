%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun,ats,bts]=IBPL(var,ngmar,aa,maxiteropt,at,bt,obj,stopindex,r,tk,btmax,rho)
%% initialization algorithm
loss=[];
timerun=[0];
tks=1;
num=length(var);
N=size(var{1},2);
for i=1:num
    varze{i}=zeros(size(ngmar,i),N);
end
gamma=1.01;
if(obj==2)
    tks=(1+sqrt(1+4*tks^2))/2;
    bt=(tks-1)/(tks);
end
varK=var;





at=min(gamma*bt,btmax);
bts=[];
ats=[];


returnloss=norm(ngmar)^2;
loss(1)=compute(var,ngmar);
t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
for j=1:num
    varz{j}=var{j}+at*(var{j}-varK{j});
    vv{j}=var{j}+bt*(var{j}-varK{j});
end
varK=var;
var=IBPLupdate(vv,varz,varK,varze,num,ngmar,aa,r);
loss(i+1)=compute(var,ngmar);


%% Judging whether to extrapolate
checkside=IBPLjudge(vv,varz,var,loss(i+1),loss(i),rho);
if(checkside==1)
    var=varK;
    var=IBPLupdate(var,var,varK,varze,num,ngmar,aa,r);
    loss(i+1)=compute(var,ngmar);
end

if(obj==1 && checkside==1)
    bt=bt/tk; 
elseif (obj==1 && checkside==0)
    bt=min(bt*tk,btmax);
end


if(obj==2)
    tks=(1+sqrt(1+4*tks^2))/2;
    bt=(tks-1)/(tks);
end

at=min(gamma*bt,btmax);
bts(i)=bt;
ats(i)=at;

%% Check if termination condition is met
fprintf("IBPL\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
t2=clock;
timerun(i+1)=etime(t2,t1);
ReeK=loss(i)/returnloss;
Ree=loss(i+1)/returnloss;
Res=abs(Ree-ReeK);
fprintf("Rel：%d\n",Res);
stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end
end
end


function checkside=IBPLjudge(vv,varz,var,lossplus,loss,rho)
sumnorm=0;
checkside=0;
for i=1:length(var)
    sumnorm=sumnorm+norm(vv{i}-var{i},'fro')^2+norm(varz{i}-var{i},'fro')^2;
end
if(lossplus>loss-rho/2*sumnorm)
      checkside=1;
end
end












