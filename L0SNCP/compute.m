%% Calculate the value of the objective function
function loss=compute(var,ngmar)
X=ktensor(var);
S=tensor(X);
loss=norm(ngmar-S)^2/2;
end

