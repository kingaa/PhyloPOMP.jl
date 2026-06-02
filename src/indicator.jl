"""
    @indicator condition value

Returns `value` if `condition=true`, `zero(Float64)` otherwise.
"""
macro indicator(condition,value)
    expr = quote
        thunk = () -> $value
        if $condition
            thunk()
        else
            zero(Float64)
        end
    end
    esc(expr)
end
