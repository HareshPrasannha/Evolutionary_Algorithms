function y=Mutate(x)

    n=numel(x);
    
    i=randsample(n,2);		%% Identifying to random index in the particle to mutate
    i1=i(1);
    i2=i(2);

    y=x;
    y([i1 i2])=x([i2 i1]);		%% Perform mutation on the particle

end
