%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
%%  This algorithm comes from the paper: "Proximal gradient method with extrapolation and line search for a class of non-convex and non-smooth problems"
function [var,loss,timerun,bts]=PGels(var,ngmar,aa,maxiteropt,stopindex,r)
%% initialization algorithm

C=10^-4;
eta=0.8;
Ntrace=2;
mumin=10^-6;
delta=0.1;

loss=[];
timerun=[0];
num=length(var);
L=ones(1,num);
N=size(var{1},2);

for i=1:num
    varze{i}=zeros(size(ngmar,i),N);
end

tk=1;
bts=[];
varK=var;

loss(1)=compute(var,ngmar);
lossH(1)=loss(1);
returnloss=norm(ngmar)^2;


wk=zeros(1,num);
for j=1:num
wk(j)=(tk-1)/(tk);
end
t1=clock;
for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
vv=var;
% mumaxK=mumax;
mumax=(L+2*C)/(1-delta);

for j=1:num
    %% This is the extrapolation parameter update strategy adopted by the BPL algorithm.
    wk(j)=min(wk(j),sqrt(delta*(mumax(j)-L(j))*mumax(j)/(4*(mumax(j)+L(j))^2)));
    %% Update parameters
    vv{j}=var{j}+wk(j)*(var{j}-varK{j});
    varK{j}=var{j};
    [V,L(j)]=grad(vv,ngmar,j,num,r);
    var{j}=PROX(varze{j},V,aa(j));
    vv{j}=var{j};
end

loss(i+1)=compute(var,ngmar);
lossH(i+1)=loss(i+1);
sumnorm=0;
for j=1:length(var)
    temp=norm(var{j}-varK{j},'fro')^2;
    sumnorm=sumnorm+temp;
    lossH(i+1)=lossH(i+1)+delta*mumax(j)*temp/4;
end



%% Judging whether to extrapolate
if(lossH(i+1)>lossH(i)-C/2*sumnorm)
    var=varK;
    for j=1:num
    [V,L(j)]=grad(var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},V,aa(j));
    end
    loss(i+1)=compute(var,ngmar);
    lossH(i+1)=loss(i+1);
    for j=1:length(var)
    temp=norm(var{j}-varK{j},'fro')^2;
    lossH(i+1)=lossH(i+1)+delta*mumax(j)*temp/4;
    end
end



%% Check if termination condition is met
fprintf("PGels\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));
end
bts{i}=wk;
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
tk=(1+sqrt(1+4*tk^2))/2;
for j=1:num
wk(j)=(tk-1)/(tk);
end
end

end












