
module DynamicIntervals

export DInt, +, *, deriv, duole_theorem, range_contraction

using LinearAlgebra
using Statistics  # 用于数值导数

# 动态区间数类型（论文第一部分定义1）
struct DInt{T<:Real}
    b::Vector{T}       # 中心值 b(t)
    δm::Vector{T}      # 下界偏差 δ⁻(t)（非正）
    δp::Vector{T}      # 上界偏差 δ⁺(t)（非负）
    ε::Dict{Symbol,Vector{T}}  # 依赖跟踪：键为变量名，值为系数αⱼ(t)
end

# 构造函数：确保δm≤0，δp≥0，且ε中向量长度与b一致
function DInt(b::Vector{T}, δm::Vector{T}, δp::Vector{T}, ε::Dict{Symbol,Vector{T}}=Dict()) where T<:Real
    @assert length(b) == length(δm) == length(δp) "b、δm、δp维度必须一致"
    @assert all(δm .≤ 0) "δm必须为非正值"
    @assert all(δp .≥ 0) "δp必须为非负值"
    # 检查ε中所有向量长度与b一致
    for (k, v) in ε
        @assert length(v) == length(b) "ε[$k]长度与b不一致"
    end
    DInt{T}(b, δm, δp, ε)
end

# 标量构造（简化1维情况）
DInt(b::T, δm::T, δp::T, ε::Dict{Symbol,Vector{T}}=Dict()) where T<:Real = 
    DInt([b], [δm], [δp], Dict(k => [v[1]] for (k, v) in ε))

# 基础加法（论文第一部分运算规则）
function Base.:+(a::DInt{T}, b::DInt{T}) where T
    @assert length(a.b) == length(b.b) "加法：a和b维度必须一致"
    new_b = a.b + b.b
    new_δm = a.δm + b.δm
    new_δp = a.δp + b.δp
    # 合并依赖项（检查同键向量长度一致）
    new_ε = Dict{Symbol,Vector{T}}()
    for k in union(keys(a.ε), keys(b.ε))
        v1 = get(a.ε, k, zeros(T, length(a.b)))
        v2 = get(b.ε, k, zeros(T, length(b.b)))
        @assert length(v1) == length(v2) "加法：ε[$k]长度不一致"
        new_ε[k] = v1 + v2
    end
    DInt(new_b, new_δm, new_δp, new_ε)
end

# 基础乘法（一阶区间扩张+高阶包络，简化版）
function Base.:*(a::DInt{T}, b::DInt{T}) where T
    @assert length(a.b) == length(b.b) "乘法：a和b维度必须一致"
    new_b = a.b .* b.b
    # 一阶偏差：|b₁|δ₂ + |b₂|δ₁（上下界分别计算）
    new_δm = abs.(a.b) .* b.δm + abs.(b.b) .* a.δm
    new_δp = abs.(a.b) .* b.δp + abs.(b.b) .* a.δp
    # 合并依赖项（检查同键向量长度一致）
    new_ε = Dict{Symbol,Vector{T}}()
    for k in union(keys(a.ε), keys(b.ε))
        v1 = get(a.ε, k, zeros(T, length(a.b)))
        v2 = get(b.ε, k, zeros(T, length(b.b)))
        @assert length(v1) == length(v2) "乘法：ε[$k]长度不一致"
        new_ε[k] = v1 + v2  # 简化：实际可能需要更复杂的乘法规则
    end
    DInt(new_b, new_δm, new_δp, new_ε)
end

# 导数运算（论文第一部分规则）
function deriv(n::DInt{T}, t::Vector{T}) where T  # t为时间点
    @assert length(t) == length(n.b) "时间点与数据长度必须一致"
    # 数值导数（处理长度为1的情况）
    db = gradient(n.b, t)
    dδm = gradient(n.δm, t)
    dδp = gradient(n.δp, t)
    # 依赖项导数（确保长度一致）
    new_ε = Dict(k => gradient(v, t) for (k, v) in n.ε)
    DInt(db, dδm, dδp, new_ε)
end

# 辅助函数：计算数值梯度（修复长度为1的情况）
function gradient(y::Vector{T}, x::Vector{T}) where T
    n = length(y)
    @assert n == length(x) "y和x长度必须一致"
    n == 0 && error("输入向量不能为空")
    dy = similar(y)
    if n == 1
        # 长度为1时导数定义为0（或根据实际需求修改）
        dy[1] = zero(T)
    else
        for i in 1:n
            if i == 1
                dy[i] = (y[2] - y[1]) / (x[2] - x[1])  # 前向差分
            elseif i == n
                dy[i] = (y[n] - y[n-1]) / (x[n] - x[n-1])  # 后向差分
            else
                dy[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])  # 中心差分
            end
        end
    end
    return dy
end

# 多乐定理（紧致误差界，论文第一部分定理4）
function duole_theorem(n::DInt{T}, L::T, t::T, n0::DInt{T}, t0::T) where T
    @assert length(n0.b) == 1 "初始状态n0应为标量（1维）"  # 简化假设
    @assert t > t0 "时间t必须大于t0"
    haus0 = max(norm(n0.δm), norm(n0.δp))  # 初始Hausdorff距离
    # 积分项：∫₀ᵗ e^(L(t-τ))·max(δ⁻,δ⁺)dτ（确保积分区间长度≥2）
    τ = range(t0, t, length=max(2, 100))  # 至少2个点避免空向量
    max_δ = [max(norm(n.δm), norm(n.δp)) for _ in τ]
    integral = trapz(τ, exp.(L .* (t .- τ)) .* max_δ)
    bound = exp(L * (t - t0)) * haus0 + integral  # 修正：用t-t0而非t
    return bound
end

# 梯形法积分（辅助函数，确保输入长度≥2）
function trapz(x::Vector{T}, y::Vector{T}) where T
    @assert length(x) == length(y) ≥ 2 "x和y必须至少有2个元素"
    0.5 * dot((x[2:end] - x[1:end-1]), (y[1:end-1] + y[2:end]))
end

# 范围收缩验证（论文第一部分定理2）
function range_contraction(n::DInt{T}, δ_threshold::T) where T
    δ_norm = max(norm(n.δm), norm(n.δp))
    is_contracted = δ_norm < δ_threshold
    error_order = δ_norm  # O(‖δ‖)
    return (is_contracted, error_order)
end

end # module
