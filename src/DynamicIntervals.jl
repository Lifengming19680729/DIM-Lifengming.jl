module DynamicIntervals

export DInt, +, -, *, deriv, duole_theorem, range_contraction, lower_bound, upper_bound

using LinearAlgebra

"""
    DInt(b, δm, δp)
动态区间数类型。
- `b`: 中心值（标量或向量）
- `δm`: 下界偏差（非正）
- `δp`: 上界偏差（非负）
表示区间 [b + δm, b + δp]。
"""
struct DInt{T<:Real}
    b::T
    δm::T
    δp::T
    function DInt(b::T, δm::T, δp::T) where T<:Real
        @assert δm <= 0 "δm must be non-positive"
        @assert δp >= 0 "δp must be non-negative"
        new{T}(b, δm, δp)
    end
end

# 便捷构造函数：对称偏差
DInt(b, δ) = DInt(b, -δ, δ)

lower_bound(x::DInt) = x.b + x.δm
upper_bound(x::DInt) = x.b + x.δp

# 加法规则: (b1; δm1, δp1) + (b2; δm2, δp2) = (b1+b2; δm1+δm2, δp1+δp2)
function Base.:+(a::DInt, b::DInt)
    return DInt(a.b + b.b, a.δm + b.δm, a.δp + b.δp)
end

# 减法规则: (b1; δm1, δp1) - (b2; δm2, δp2) = (b1-b2; δm1-δp2, δp1-δm2)
function Base.:-(a::DInt, b::DInt)
    return DInt(a.b - b.b, a.δm - b.δp, a.δp - b.δm)
end

# 乘法规则 (一阶近似): (b1; δm1, δp1) * (b2; δm2, δp2) ≈ (b1*b2; b1*δm2 + b2*δm1, b1*δp2 + b2*δp1)
function Base.:*(a::DInt, b::DInt)
    new_b = a.b * b.b
    new_δm = a.b * b.δm + b.b * a.δm
    new_δp = a.b * b.δp + b.b * a.δp
    return DInt(new_b, new_δm, new_δp)
end

# 导数运算（数值微分，中心差分）
function deriv(n::Function, t::Real; h=1e-5)
    n_plus = n(t + h)
    n_minus = n(t - h)
    db = (n_plus.b - n_minus.b) / (2h)
    dδm = (n_plus.δm - n_minus.δm) / (2h)
    dδp = (n_plus.δp - n_minus.δp) / (2h)
    return DInt(db, dδm, dδp)
end

# 多乐定理（紧致误差界）
function duole_theorem(n::DInt, L::Real, t::Real, n0::DInt, t0::Real)
    haus0 = max(abs(n0.δm), abs(n0.δp))
    integral_term = (haus0 / L) * (exp(L * (t - t0)) - 1)
    bound = exp(-L * (t - t0)) * haus0 + integral_term
    return bound
end

# 范围收缩验证
function range_contraction(n::DInt, δ_threshold::Real)
    δ_max = max(abs(n.δm), abs(n.δp))
    return (δ_max < δ_threshold, δ_max)
end

end # module
