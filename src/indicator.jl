"""
    @indicator condition value type

Returns `value` if `condition=true`, `zero(type)` otherwise.
By default, `type = Float64`.
"""
macro indicator(condition, value, type = Float64)
    esc(:(($condition) ? $value : zero($type)))
end
