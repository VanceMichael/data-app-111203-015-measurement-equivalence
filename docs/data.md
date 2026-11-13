# 数据契约

尺寸检测台账保存新旧工艺配对测量，按产线提供差值均值和有效配对率。

主测量 `new_value` 的单位为微米，主体字段是 `sample`，核算维度是 `line`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `observation_id` | str | 必填 |
| `sample` | str | 必填 |
| `line` | str | 必填 |
| `day` | int | 必填 |
| `new_value` | float | 必填 |
| `old_value` | float | 允许空值，空值不是零 |
| `operator` | str | 必填 |
| `lot` | str | 必填 |
| `valid` | int | 必填 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。

## 业务字段语义

new_value 与 old_value 均为微米；old_value 可空。valid=1 表示可用配对，valid=0 表示不可用于配对比较，两者都不进入配对推断。sample 可跨操作员和日期重复测量，lot 是样品批次，line 为产线。

## 等效性研究中的使用

- 配对差值取 `new_value - old_value`，仅用于 `valid=1` 且旧值非空的观测。
- `sample` 是重复测量的分析单位：同一 sample 的多条观测先合并为差值均值，
  再做配对 TOST，重复送样不按独立样本计数。
- `lot` 与 `operator` 作为聚类维度处理相关性；`sample` 必须各自只属于一个
  `lot`。跨批的 sample 无法定义批次聚类，该判定记为不可估计；跨线的
  sample 会记录在 `scope_issues` 并阻止等效声明，但不进入任何产线亚组。
- `line` 是预设亚组，计划中列出的每条产线都要判定；不能事后挑选通过的产线。

字段类型与空值规则仍遵循上表，额外列不进入当前核算。
