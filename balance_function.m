function n = balance_function(a,dhdx,g,rho,x,A,H_obs)

    fun = @(n) sqrt(sum((((-n+2)./(2*A).*abs(dhdx).^(1-n).*dhdx.^(-1).*a).^1/(n+2)- H_obs).^2));
    n0 = 4;
    fun(1:4)
    plot(1:4, fun(1:4))
    [n] = fminunc(fun,n0);

end