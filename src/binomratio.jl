macro binomratio(n, ℓ, r::Integer, s::Integer, type::Type = Float64)
    expr = if 0 ≤ s ≤ r
        cond = :($s ≤ $ℓ ≤ $n && $r ≤ $n)
        if r == 0
            :(@indicator($cond, one($type), $type))
        elseif r ==1 && s == 0
            :(@indicator($cond, ($n-$ℓ)/$n, $type))
        elseif r == 1 && s == 1
            :(@indicator($cond, inv($n), $type))
        elseif r == 2 && s == 0
            :(@indicator($cond, ($n-$ℓ)/($n)*($n-$ℓ-1)/($n-1), $type))
        elseif r == 2 && s == 1
            :(@indicator($cond, 2*($n-$ℓ)/($n)/($n-1), $type))
        elseif r == 2 && s == 2
            :(@indicator($cond, 2/($n)/($n-1), $type))
        else
            :(@indicator($cond, binomial($n-$ℓ,$(r-s))/binomial($n,$r), $type))
        end
    else
        :(zero($type))
    end
    esc(expr)
end

macro binomratio(n, ℓ, r, s, type::Type = Float64)
    if n.head == ℓ.head == r.head == s.head == :tuple &&
        length(n.args) == length(ℓ.args) == length(r.args) == length(s.args)
        h = map(n.args, ℓ.args, r.args, s.args) do a,b,c,d
            :(@binomratio($a,$b,$c,$d))
        end
        esc(Expr(:call, :*, h...))
    else
        error("invalid @binomratio construction")
    end
end
