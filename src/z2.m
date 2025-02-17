function out = z2(params)


t = params.t;

if ( t < 0.5 )
    out = (2*t);
else
    out = t+0.5;
end



end