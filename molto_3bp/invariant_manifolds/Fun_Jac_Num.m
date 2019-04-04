function Jnum = Fun_Jac_Num(funcion,t,Y,dh)
    
Jnum = zeros(length(Y),length(Y));

  
    for cont = 1:1:length(Y)

        YM       = Y;
        YM(cont) = YM(cont) + dh/2;
        fM       = feval(funcion,t,YM);

        Ym       = Y;
        Ym(cont) = Ym(cont) - dh/2;
        fm       = feval(funcion,t,Ym);     

        Jnum(:,cont) = (fM - fm)/dh;        
    end
end