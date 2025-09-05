function lossrecord=valueplot(datasss,timeend,trigger)
warning("off")
trinum=trigger;
trihigh=[7,8,9];
% trihigh=[8,9];

iter=0.1;
x1 = 0:iter:timeend;
percent=length(x1)/2;
% percent=length(x1)/2.9;
% percent=length(x1)/3;
% percent=length(x1)/2;

ss=20;


if(isempty(timeend))
    timeend=400;
end
color=["-o","-p","-d","-x","->","-*","-+","-^","-<"];
RGB={[0 0.4470 0.7410],[0 0 0.8],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[1 0 0],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0.8500 0.3250 0.0980]};








lossdatas=[];
tr=[];
siters=[];
lsdata=[];
for k=1:length(datasss)

datass=datasss{k};

for i=1:length(datass)-1
    datas=datass{i};

    lsdata{i}=datas{1};
    rtdata{i}=datas{2};
end
lossdata{k}=lsdata;
timedata{k}=rtdata;
datas=[];
lsdata=[];
rtdata=[];
end




for k=1:length(lossdata)
datas=lossdata{k};
tdatas=timedata{k};




for j=1:length(datas)
datatemp=datas{j};
tdatatemp=tdatas{j};
index=find(tdatatemp>timeend);
inci=min(index);
data=datatemp(1:inci);
tdata=tdatatemp(1:inci);
if((tdata(end-1)-timeend)<1e-5)
    tdata(end)=timeend+0.01;
else
    tdata(end)=timeend;
end
y1 = interp1(tdata,data,x1,'linear');
datainter{j}=y1;
end
datainters{k}=datainter;
datainter=[];
data=[];
tdata=[];
end




a=datainters{1};
for i=1:length(a)
    for k=1:length(datainters)
    datainter=datainters{k};
    data=datainter{i};
    lossrecord_global(k,i)=data(end);
    end
end






loss=[];
a=datainters{1};
for i=1:length(a)
    loss=zeros(1,size(0:iter:timeend,2));
    for k=1:length(datainters)
    datainter=datainters{k};
    data=datainter{i};
    lossrecord(k,i)=data(end);
    loss=loss+data;
    end
    loss=loss/length(datainters);
    losssum{i}=loss;
end






figure(1)
subplot(1,2,1);
mess={};
for s=1:length(trinum) 
lossdata=losssum{s};

    
    maker_idx = 1:ss:length(lossdata);
    semilogy(x1,lossdata,color(s),'linewidth',3,'MarkerIndices',maker_idx,'MarkerSize', 8,'color',RGB{s});
    xlabel('Time(s)','FontSize',30);
    ylabel('Objective funciton value','FontSize',30);
    if(trigger(s)==1)
    mes='iPALM';
    elseif(trigger(s)==2)
    mes='BPL';
    elseif(trigger(s)==3)
     mes='IBPG';
    elseif(trigger(s)==4)
    mes='TITAN';
    elseif(trigger(s)==5)
    mes='ABPL';
    elseif(trigger(s)==6)
    mes='PGels';
   elseif(trigger(s)==7)
    mes='APGL';  
    elseif(trigger(s)==8)
    mes='IBPL'; 
    elseif(trigger(s)==9)
    mes='IBPL$^{+}$';  
    end
    ls=length(mess);
    mess{ls+1}=mes;
    mes={};

    hold on;

end

prettyAxes().gbase2()
h=legend(mess,'Interpreter','latex');
set(gca,'FontSize',45);
set(h,'FontSize',28)
xlabel('Time (seconds)','FontSize',45);
ylabel('Objective funciton value','FontSize',45);


set (gca,'position',[0.108,0.178,0.445,0.7390]);
set (h,'position',[0.14,0.583,0.24,0.24]);

mess={};
subplot(1,2,2);
[~,ind]=intersect(trigger,trihigh);

for s=1:length(ind)  

    k=ind(s);
    lossdata=losssum{k};

    
    maker_idx = 1:ss:length(lossdata);

    semilogy(x1(ceil(percent):length(x1)),lossdata(ceil(percent):length(x1)),color(k),'linewidth',2.6,'MarkerIndices',maker_idx,'color',RGB{k},'MarkerSize', 8);
    if(trigger(k)==1)
    mes='iPALM';
    elseif(trigger(k)==2)
    mes='BPL';
    elseif(trigger(k)==3)
     mes='IBPG';
    elseif(trigger(k)==4)
    mes='TITAN';
    elseif(trigger(k)==5)
    mes='ABPL';
    elseif(trigger(k)==6)
    mes='PGels';
   elseif(trigger(k)==7)
    mes='APGL';  
    elseif(trigger(k)==8)
    mes='IBPL'; 
    elseif(trigger(k)==9)
    mes='IBPL$^{+}$'; 
    end
    ls=length(mess);
    mess{ls+1}=mes;
    hold on;
end
 hold on;
prettyAxes().gbase2()
% h=legend(mess,'Interpreter','latex');
set(gca,'FontSize',25);
% set(h,'FontSize',25)
set (gca,'position',[0.41,0.507,0.13,0.41]);
end






