%% Read data set function
function [ngmar,num]=readfile(i)
   if(i==1) 
   E=load('.\Data\breastmnist_224.mat');
   A1=E.train_images;
   A2=E.test_images;
   A3=E.val_images;
   ngmar= double(cat(1,A1,A2,A3));
   num=length(size(ngmar));
   elseif(i==2) 
   E=load('.\Data\PosteriorEmssion.mat');
   ngmar= E.tensor_nc_var;
   num=length(size(ngmar));

    end 
end