%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
%%  This algorithm comes from the paper: "An Accelerated Block Proximal Framework with Adaptive Momentum for Nonconvex and Nonsmooth Optimization"
function [var,loss,timerun,bts]=ABPL(var,ngmar,aa,maxiteropt,bt,obj,stopindex,r,tk,btmax)
%% initialization algorithm
loss=[];
timerun=[0];
% tks=1;
num=length(var);

for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end
bts=[];
varK=var;

loss(1)=compute(var,num,ngmar);
returnloss=norm(ngmar,"fro")^2;


t1=clock;
for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
for j=1:num
    vv{j}=var{j}+bt*(var{j}-varK{j});
end
varK=var;
var=ABPLupdate(vv,varze,num,ngmar,aa,r);
loss(i+1)=compute(var,num,ngmar);

%% Judging whether to extrapolate
[bt,var,loss(i+1)]=ABPLjudge(var,num,loss(i+1),loss(i),ngmar,varK,varze,aa,vv,bt,tk,obj,btmax,r);
bts(i)=bt;


%% Check if termination condition is met
fprintf("ABPL\n");
t2=clock;
timerun(i+1)=etime(t2,t1);
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
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




