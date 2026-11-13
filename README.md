# measurement-equivalence

尺寸检测台账保存新旧工艺配对测量，按产线提供差值均值和有效配对率。
`summary` 只做核对；工艺上线前的等效判定必须走独立的等效性研究流程：
先登记界值固化计划，再用独立数据做配对 TOST，详见
[`docs/equivalence.md`](docs/equivalence.md)。

## 使用

运行需要 Python 3.11 或更新版本，没有第三方运行依赖。数据在本地处理，不请求外部接口。

```bash
python3 -m analysis_app summary
python3 -m analysis_app summary --input data/observations.csv --output output/summary.json

python3 -m analysis_app plan --config examples/config/equivalence.json \
    --output examples/plans/equivalence.plan.json
python3 -m analysis_app check-plan --plan examples/plans/equivalence.plan.json
python3 -m analysis_app analyze --plan examples/plans/equivalence.plan.json \
    --input data/observations.csv --output examples/results/equivalence.json
python3 -m analysis_app power --plan examples/plans/pilot-power.plan.json \
    --input data/observations.csv --effect 0.17 --design-effect 4
python3 -m analysis_app examples
```

等效性结果只有 `equivalent`、`not_proven`、`unavailable` 三种；界值以原始
微米单位在分析前写入计划，计划指纹会随计划内容变化，旧分析不会被新计划
覆盖。

标准输出是完整 JSON，输出文件使用临时文件原子替换；失败退出码为 2，不覆盖既有结果。空数据集保留零计数与空分组，缺失数值不进入数值分母。数据指纹由按 observation_id 排序后的解析记录计算。

## 验证与容器

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q analysis_app
python3 -m analysis_app examples
docker build -t measurement-equivalence .
docker run --rm measurement-equivalence
```

容器内也可执行 `python -m unittest discover -s tests -v`。镜像构建阶段取得 Python 镜像，产品启动、测试和演示均无需下载。小规模样例的完整单元测试和 summary 可以连续执行，不需要等待定时任务。

## 代码与数据

- `analysis_app/io.py`：有类型的 CSV 读取、完整性校验及内容指纹。
- `analysis_app/schema.py`：字段类型、主体与聚合维度。
- `analysis_app/domain.py`：当前业务核算函数。
- `analysis_app/summary.py`：确定性的分组与数值口径。
- `analysis_app/plans.py`：等效性计划的固化、规范化和指纹校验。
- `analysis_app/equivalence.py`：配对 TOST、聚类稳健方差与 Holm 家族判定。
- `analysis_app/samplesize.py`：基于先导方差和给定偏移的功效曲线。
- `analysis_app/stats.py`：标准库实现的 t 分布、TOST 与 Holm。
- `analysis_app/examples.py`：可复算案例的登记与核对。
- `analysis_app/cli.py`：本地命令及原子输出。
- `docs/data.md`：数据契约。
- `docs/equivalence.md`：等效性研究流程与统计口径。
- `examples/`：合成配置、固化计划、案例数据与复算清单。

演示数据全部为本地合成，不对应现实机构、设备或人员；不能用来推断真实业务效果。
