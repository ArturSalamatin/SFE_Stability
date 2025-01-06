function out = dXdt(t, X, r)

out = 1./X.*(r+(1-r)*X)./(r+(1-r)*sqrt(2*t));

end