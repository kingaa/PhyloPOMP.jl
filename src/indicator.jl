"""
    @indicator condition value type

Returns `value` if `condition=true`, `zero(type)` otherwise.
By default, `type = Prob`.
"""
macro indicator(condition, value, type = Prob)
    esc(:(($condition) ? $type($value) : zero($type)))
end
