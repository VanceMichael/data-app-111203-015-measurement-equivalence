# 等效性研究流程

本流程用于判断新旧检测工艺在预先登记的微米界值内是否等效。它与
`summary` 的差值汇总相互独立：汇总只描述数据，等效结论只能由固化
计划、独立数据和完整判定链产生。

## 1. 分析前登记计划

界值以原始微米单位写入配置，不允许在分析阶段缩放或临时修改：

```json
{
  "lower": -0.5,
  "upper": 0.5,
  "alpha": 0.05,
  "lines": ["L0", "L1", "L2"],
  "cluster_dimensions": ["lot", "operator"],
  "require_multiple_lots": true
}
```

```bash
python3 -m analysis_app plan --config examples/config/equivalence.json \
    --output examples/plans/equivalence.plan.json
python3 -m analysis_app check-plan --plan examples/plans/equivalence.plan.json
```

`plan` 固化单位、界值、alpha、分析单位（sample）、聚类维度、预设产线、
Holm 家族和跨批范围要求，并写入 `plan_fingerprint`。计划文件一旦登记
不能覆盖；`plan` 对已存在路径直接失败。任何字段被改动后指纹校验都会
失败，分析拒绝运行。

## 2. 用独立数据判定

```bash
python3 -m analysis_app analyze --plan examples/plans/equivalence.plan.json \
    --input data/observations.csv --output examples/results/equivalence.json
```

- `valid=0`、旧值为空的观测不进入配对推断；`valid` 只接受 0 或 1。
- 重复送样按 `sample` 合并成分析单位，每条观测的残差影响按该 sample
  的观测数均分，重复测量不会被当成独立样本。
- `lot` 与 `operator` 的相关性用多向聚类稳健方差处理：先对 lot、
  operator 及其交集分别求和，再以交替符号合并，每个聚类使用
  `G/(G-1)` 修正，自由度取可用聚类数减一。
- 每个判定给出两个单侧检验：`p_lower` 对应 `H0: μ <= lower`，
  `p_upper` 对应 `H0: μ >= upper`。`p_value` 取两者较大值，
  `interval_1_minus_2alpha` 是与之对应的双侧 1-2α 区间。
- 总体与计划中列出的每一条产线共享一个 Holm 家族。不能先看哪条
  产线通过再决定发布哪条；只有家族内全部判定通过且范围核对满足，
  研究结论才是 `equivalent`。输出中的 `interval_1_minus_2alpha` 是
  单个判定的未校正区间，`p_adjusted` 才是家族校正后的判定依据。

## 3. 三种结果

| 判定 | 含义 |
| --- | --- |
| `equivalent` | 全部计划判定通过，区间落在界值内，且范围核对通过 |
| `not_proven` | 至少一个判定未通过，或只有一批材料等前提不满足 |
| `unavailable` | 有效配对不足、方差为零、聚类不足，不能给出判定 |

`unavailable` 不是“通过”，也不是“不等效”，原因是“分析单位不足以
估计方差”而不是数据支持了某个结论。单一批次时 lot 方差分量被明确
列入 `omitted_clusters`，即使其余检验显著，研究结论仍为 `not_proven`，
不能声称跨批适用。计划外出现的产线、计划内缺失的产线、跨批或跨线的
sample 都会记录在 `scope_issues`；跨批的 sample 使该判定不可估计，
其余范围问题阻止等效结论。

## 4. 样本量探索

```bash
python3 -m analysis_app power --plan examples/plans/pilot-power.plan.json \
    --input data/observations.csv --effect 0.17 --sizes 12,24,36,60 \
    --design-effect 4
```

`power` 用先导数据的 sample 级差值方差和用户显式给出的真实偏移
（`--effect`，微米）计算功效曲线，`--effect` 缺失时直接失败，避免用
观测到的差值冒充预期效应。`--design-effect` 可把聚类造成的方差放大
计入。输出中的 `smallest_size_at_target` 是达到 80% 功效的最小试算
规模；它只是设计输入，不能写进判定结果，也不能用事后功效替代等效
检验。

## 5. 可复算例子

```bash
python3 -m analysis_app examples
```

清单 `examples/cases.json` 登记了没有有效配对、零方差、均值正好碰界、
上下界不对称、重复送样、只有一批材料和一条产线失败七个案例。命令会
重新固化内联计划、重新分析，并逐项核对结论、估计值、自由度、校正
p 值、区间与聚类信息；任何不符都会以退出码 2 失败。

全部数据为本地合成，不能据此推断真实设备或工艺效果。
