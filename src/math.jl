#=
Acessory Math Functions
=#

"""
Pochammer symbol
...
# Arguments
- `c::Int`: c integer numbers.
- `m::Int`: m integer numbers.
...
"""
function Poch(c::Int, m::Int)
    # Pochammer Symbol
    c_i = 1
    for i = 0:m-1
        c_i = c_i * (c + i)
    end
    return c_i
end

"""
Kroenecker delta
"""
δ(x, y) = ==(x, y)
