# 动态数学体系（DynamicIntervals.jl）

实现论文《原始场（混沌场）理论》中的动态数学体系（DIM），支持时变不确定性的区间运算。

## 核心功能
- 动态区间数（`DInt`）类型，包含中心值、上下偏差及依赖跟踪
- 基础运算：加法、乘法、导数
- 多乐定理（紧致误差界）与范围收缩验证

## 安装
```julia
] add https://github.com/[你的GitHub用户名]/DynamicIntervals.jl
