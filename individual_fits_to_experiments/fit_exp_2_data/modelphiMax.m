function output= modelphiMax(xEst,thours,NG0)
    b0=NG0;%initial antibiotic concentration which is chosen as the time 0 value
    options = odeset('AbsTol',1e-13,'RelTol',1e-11);
   
    [~,values] = ode45(@Equations,thours,b0,options);
 
    function s=Equations(~,y)
     
    phi_max=xEst;%replication rate  of unattached bacteria
   
    s=zeros(1,1);
    s(1)=phi_max*y(1);    
   
    end
 
output=values(:,1);

 end