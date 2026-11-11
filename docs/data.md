# 数据契约

尺寸检测台账保存新旧工艺配对测量，按产线提供差值均值和有效配对率。

主测量 `new_value` 的单位为微米，主体字段是 `sample`，核算维度是 `line`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `observation_id` | str | 必填 |
| `sample` | str | 必填 |
| `line` | str | 必填 |
| `day` | int | 允许空值，空值不是零 |
| `new_value` | float | 必填 |
| `old_value` | float | 允许空值，空值不是零 |
| `operator` | str | 允许空值，空值不是零 |
| `lot` | str | 允许空值，空值不是零 |
| `valid` | int | 允许空值，空值不是零 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。
