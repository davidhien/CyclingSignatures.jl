# custom finite field type
function is_prime(n::Int)
    if iseven(n) || n < 2
        return n == 2
    else
        p = 3
        q = n ÷ p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n ÷ p
        end
        return true
    end
end
# The following causes an error if you try to use anything other than Int as a modulus.
is_prime(::Any) = false

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where {M}
    is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
    i = i % M
    return i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) = i & 1

"""
    FF{M} <: Integer

`FF{M}` is the default field used in CyclingSignatures. It is a representation of a finite field
``\\mathbb{Z}_M``, integers modulo small, prime `M`. Supports field arithmetic and can be
converted to integer with `Int`.


# Example

```jldoctest
julia> FF{3}(5)
2 mod 3

julia> FF{3}(5) + 1
0 mod 3

```
"""
struct FF{M} <: Integer
    value::Int

    # Check mod allows construction when you know you don't need to mod the number.
    function FF{M}(value::Integer, check_mod=true) where {M}
        if check_mod
            return new{M}(mod_prime(value, Val(M)))
        else
            return new{M}(value)
        end
    end
end
FF{M}(i::FF{M}) where {M} = i

Base.Int(i::FF) = i.value

Base.show(io::IO, i::FF{M}) where {M} = print(io, Int(i), " mod ", M)

for op in (:+, :-, :*)
    @eval (Base.$op)(i::FF{M}, j::FF{M}) where {M} = FF{M}($op(Int(i), Int(j)))
end

Base.:/(i::FF{M}, j::FF{M}) where {M} = i * inv(j)
Base.:-(i::FF{M}) where {M} = FF{M}(M - Int(i), false)
Base.zero(::Type{FF{M}}) where {M} = FF{M}(0, false)
Base.one(::Type{FF{M}}) where {M} = FF{M}(1, false)
Base.sign(i::M) where {M<:FF} = ifelse(iszero(i), zero(M), one(M))

Base.promote_rule(::Type{FF{M}}, ::Type{<:FF}) where {M} = Union{}
function Base.promote_rule(::Type{FF{M}}, ::Type{I}) where {M,I<:Integer}
    if Base.promote_type(I, Int128) === Int128 || Base.promote_type(I, Int128) === UInt128
        return FF{M}
    else
        return Union{}
    end
end

Base.inv(i::FF{M}) where {M} = FF{M}(invmod(Int(i), M), false)