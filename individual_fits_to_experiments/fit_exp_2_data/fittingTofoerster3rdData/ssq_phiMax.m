function SSE=ssq_phiMax(parameters,x)
    time=x(:,1);
    data1=x(:,2);
        
    NG0=data1(1); 
   
    fitted=modelphiMax(parameters,time,NG0);
   
    error1=log10(fitted(:,1))-log10(data1);
    SSE=sum(error1.^2);
        
end