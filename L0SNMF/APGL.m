%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
%%  This algorithm comes from the paper: "Sparse nonnegative CP decomposition with graph regularization and $\ell_{0}$-constraints for clustering"
function [var,loss,timerun,bts]=APGL(var,ngmar,aa,maxiteropt,stopindex,btmax,rho)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));
for i=1:num
    varze{i}=zeros(size(var{i}));
end
r=1.5;
btmax=0.99;
LK=zeros(1,num);
L=ones(1,num);
tk=1;
bts=[];
wk=zeros(1,num);
varK=var;
loss(1)=compute(var,num,ngmar);
returnloss=norm(ngmar,"fro")^2;
t1=clock;


for i=1:maxiteropt
fprintf("%d\n",i);
vv=var;

%% Update parameters
for j=1:num
    wk(j)=min(wk(j),btmax);
    vv{j}=var{j}+wk(j)*(var{j}-varK{j});
    varK{j}=var{j};
    var{j}=vv{j};
    LK(j)=L(j);
    [V,L(j)]=grad(var,ngmar,j,num,r);
    var{j}=PROX(varze{j},V,aa(j));
end
loss(i+1)=compute(var,num,ngmar);

sumnorm=0;
for j=1:length(var)
    sumnorm=sumnorm+norm(var{j}-varK{j},'fro')^2;
end

%% Judging whether to extrapolate
if(loss(i+1)>loss(i)-rho/2*sumnorm)
    var=varK;
    for j=1:num
    [V,L(j)]=grad(var,ngmar,j,num,r);
    var{j}=PROX(varze{j},V,aa(j));
    end
    loss(i+1)=compute(var,num,ngmar);
end


bts{i}=wk;
t2=clock;
timerun(i+1)=etime(t2,t1);
fprintf("APGL\n");
%% Check if termination condition is met
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
ReeK=loss(i)/returnloss;
Ree=loss(i+1)/returnloss;
Res=abs(Ree-ReeK);
fprintf("cri：%d\n",Res);
stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end
    tktemp=tk;
    tk=(1+sqrt(1+4*tk^2))/2;
    for j=1:num
        wk(j)=(tktemp-1)/(tk);
    end
end

end


