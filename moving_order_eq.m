function dydt = moving_order_eq(t,y,nu,V,eta)


dydt = [y(2);nu^(1/2)/V*y(1)^(1/2)*(eta*y(1)-y(2))^(3/2)];


end

    
    