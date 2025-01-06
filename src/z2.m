function out = z2(t, r)
out = 1/(1-r)*(sqrt(2*t) - r/(1-r)*log(1+(1-r)/r*sqrt(2*t)));
end
