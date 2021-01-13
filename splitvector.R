

### splitvector function ###


splitvector = function(x,n)
{
    split(x,ceiling(seq_along(x)/(length(x)/n)))
}



