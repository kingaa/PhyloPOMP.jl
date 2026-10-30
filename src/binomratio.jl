macro binomratio(n, ell, r::Integer, s::Integer, type::Type = Float64)
    expr = if 0 ≤ s ≤ r
        cond = :($s ≤ $ell ≤ $n && $r ≤ $n)
        if r == 0
            :(@indicator($cond, one($type), $type))
        elseif r ==1 && s == 0
            :(@indicator($cond, ($n-$ell)/$n, $type))
        elseif r == 1 && s == 1
            :(@indicator($cond, inv($n), $type))
        elseif r == 2 && s == 0
            :(@indicator($cond, ($n-$ell)/($n)*($n-$ell-1)/($n-1), $type))
        elseif r == 2 && s == 1
            :(@indicator($cond, 2*($n-$ell)/($n)/($n-1), $type))
        elseif r == 2 && s == 2
            :(@indicator($cond, 2/($n)/($n-1), $type))
        else
            :(@indicator($cond, binomial($n-$ell,$(r-s))/binomial($n,$r), $type))
        end
    else
        :(zero($type))
    end
    esc(expr)
end

macro binomratio(n, ell, r, s, type::Type = Float64)
    if n.head==:tuple && ell.head==:tuple &&
        r.head==:tuple && s.head==:tuple &&
        length(n.args)==length(ell.args)==length(r.args)==length(s.args)
        n = n.args
        ell = ell.args
        r = r.args
        s = s.args
        h = map(n,ell,r,s) do a,b,c,d
            :(@binomratio($a,$b,$c,$d))
        end
        esc(Expr(:call,:*,h...))
    else
        error("invalid @binomratio construction")
    end
end
