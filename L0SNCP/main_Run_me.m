clearvars
clc


%% Parameter.
%   N         : Rank of Tensor CP Decomposition
%   index     : The dataset to be used, when index=1, use breastmnist_224 dataset,
%               when index=2, use PosteriorEmssion dataset
%   at,bt     : Initial extrapolation parameter values of extrapolation points of IBPL+
%   btmax     : Maximum extrapolation parameter values of IBPL+
%   t         : Decay and increment rates for extrapolation parameters for IBPL+
%   r         : Step factor
%   rho       : Parameter for the restart step
%   maxiteropt: Maximum iteration alloted to the algorithm
%   outer     : The number of times the loop runs independently
%   trigger   : Whether to enable the indicator array of each algorithm, where
%               when 1∈trigger, enable the iPALM algorithm
%               when 2∈trigger, enable the BPL algorithm
%               when 3∈trigger, enable the IBPG algorithm
%               when 4∈trigger, enable the TITAN algorithm
%               when 5∈trigger, enable the ABPL algorithm
%               when 6∈trigger, enable the PGels algorithm
%               when 7∈trigger, enable the APGL algorithm
%               when 8∈trigger, enable the IBPL+ with the non-adaptive momentum
%               when 9∈trigger, enable the IBPL+ with the adaptive momentum
%   percent   ：The proportion of non-zero elements allowed in each decomposition matrix
%   stopindex : The indicator of the stop condition.  
%               To set the specific termination condition, see the 'stopcheck' function for details.  
%               The default termination condition is: each algorithm runs for 30 seconds
%% Display
%   nonzero   ：The number of non-zero elements in each component.
%   Rel       : The difference in the objective function value between two iterations.


%% Parameter settings
warning('off');
rng('shuffle')
R=50;
index=2;
r=1.01;
bt=0.2;
at=1.01*bt;
t=1.3;
rho=10^-5;
btmax=0.9999;
outer=20;
maxiteropt=6000000;
trigger=[1,2,3,4,5,6,7,8,9];
percent=0.3;
stopindex=1;





%% Algorithm iteration starts

[ngmar,num]=readfile(index);
ngmar=tensor(ngmar);




N=R;
for i=1:length(size(ngmar))
    pere(i)=N*size(ngmar,i)*percent;
end
aa=pere;



for j=1:outer

for i=1:num
    den=size(ngmar,i)*N;
    var{i}=sprand(size(ngmar,i),N,aa(i)/(2*den));
    var{i}=full(var{i});
end


for i=1:length(trigger)
[datas{i},vars{i}]=ALGOchoose(var,ngmar,aa,maxiteropt,at,bt,trigger(i),stopindex,r,t,btmax,rho);
end
datas{length(trigger)+1}=var;
datass{j}=datas;
end





%% Output and Drawing
lossrecord=valueplot(datass,40,trigger);
Obj=(sum(lossrecord)/length(lossrecord))';
ObjErr=std(lossrecord,0,1)';
temp=sqrt(lossrecord*2)/norm(tensor(ngmar));
Rel=(sum(temp)/length(temp))';
RelErr=std(temp,0,1)';
[~,ranking]=min(lossrecord,[],2);
ranking=tabulate(ranking);






