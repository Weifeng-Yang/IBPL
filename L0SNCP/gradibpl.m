%% Gradient calculation function for IBPL
function U=gradibpl(var,varz,varK,ngmar,n,num,r)
  [Xtemp,temp]=krob2(varK,n,num,ngmar);
   a=0;
   ck=norm(temp,'fro');
   tao=1/(r*ck);
   if(temp==0) 
       while(a==0)
           a=rand(1);
       end
       tao=1/a;
   end
   mar=var{n}*temp-Xtemp;
   U=varz{n}-tao*mar;
end