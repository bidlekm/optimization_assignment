function n = balance_function(a,dhdx,g,rho,x,A,H_obs)
    fun = @(n)norm(nthroot((-n+2)./(2.*A).*abs(dhdx).^(1-n).*dhdx.^-1.*a ,n+2)- H_obs,2);
    n0 = 3;

    [n] = fminunc(fun,n0);
end