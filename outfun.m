function stop = outfun(x,optimValues,state)
global Y;
stop = false;
i=optimValues.iteration;
Y(1,i+1)=i;
Y(2:3,i+1)=x;
Y(4,i+1)=optimValues.fval;
end

