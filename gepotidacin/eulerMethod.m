function [t,solutions,te]=eulerMethod(func,t0,y0,stepSize,time,paras)
numberOfSteps=round((time-t0)/stepSize);
t=NaN(numberOfSteps,1);
solutions=NaN(numberOfSteps,6);
for j=1:numberOfSteps
    totalNG=y0(1)+y0(2)+y0(3)+y0(5);
    if totalNG <10%similar to the event function
        break
    else
        m=func(t0,y0,paras);
        y1=y0+stepSize.*m';
        t1=t0+stepSize;
        t(j)=t1;
        solutions(j,:)=y1;
        t0=t1;
        y0=y1;
    end
end
t(any(isnan(t),2), :) = [];
solutions(any(isnan(solutions),2), :) = [];
te=t(end);

end