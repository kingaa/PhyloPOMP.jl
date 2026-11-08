switchblock1(expr) = expr

switchblock1((cond, action), expr...) = Expr(:elseif, cond, action, switchblock1(expr...))

switchblock(expr...) = begin
    e = switchblock1(expr...)
    e.head = :if
    e
end
