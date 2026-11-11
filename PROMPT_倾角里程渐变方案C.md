# 方案 C：倾角沿里程连续渐变（新增第三档，原两档不动）

日期：2026-11-10
状态：**方案评审稿，未改代码**
触发问题：大平面（1030×487 这类）上所有点共用一个倾角，远端被近端绑住，卡 `JOINT_LIMIT`（腕部拧死）。

> **评审先读第 3 节。** 那节排除了一个看起来最省事、实际无效的路子（绕枪轴滚转），
> 并给出「哪些量真的能改 θ5」的推导。跳过它容易在实现阶段做无用功。

---

## 0. 一句话

新增 `RS_POSE_TILT_RAMP`（bit2）：**每行两端各自求出「满足关节余量的最小倾角」，
中点按里程线性插值**。面级常量（默认）与逐点连续（bit1）一行不改。

---

## 1. 为什么现有两个入口在大平面上都等于"没变"

| 入口 | 驱动量 | 大平面上的结果 |
|---|---|---|
| 默认：面级常量 | 面组法向均值 | 全平面一个 `(tilt, az)` |
| bit1 `RS_POSE_TILT_PERP` | `tiltNeedDeg(normal(i))` | 法向恒定 ⟹ **与默认逐位相同** |
| bit2 `RS_POSE_TILT_RAMP`（本方案） | 里程 | 两端不同 ⟹ 才有梯度 |

`tiltNeedDeg`（`robot_safety.cpp:1370`）只吃 `dot(n, tiltRefW)`。大平面法向是常量，
所以 bit1 在这类件上产生不了任何变化 —— 这是"开了逐点连续却还是拧死"的原因。

**面组判据（`robot_safety.cpp:1971`）**：相邻行法向 ≤2° 且共面 ≤2mm。
大平面全部成立 ⟹ 恒为 1 个面组 ⟹ 整面锁死一个姿态。

---

## 2. 关键约束：喷角 ≡ 倾角

`dirOf`（`robot_safety.cpp:1321`）里 `d` 是 `n` 绕切向轴转 `tiltDeg`，
所以 **枪轴与法向的夹角恰好等于倾角**：

```
sprayAngle == tiltDeg
```

后果：**加大倾角就是加大喷角**，直接吃 `maxSprayAngleDeg`（默认 30°）预算。
所以本方案的倾角必须**封顶在工艺喷角**，超了就如实报"这一端做不到"，不许偷着给 35°。

---

## 3. 杠杆分析：J5 卡死时到底哪个量有用

先排除一个想当然的路子，避免白做。

`kJoints`（`robot_safety.cpp:318`）：**J4 ±270° / J5 ±120° / J6 ±360°**，
腕部是 **Z-Y-Z 球腕**（J4 绕前臂滚转、J5 俯仰、J6 绕枪轴滚转）。

由此可得：

```
d = Rz(θ4)·Ry(θ5)·ẑ    ⟹    θ5 = arccos(d_z)
```

**θ5 只由枪轴方向 `d` 决定，与绕枪轴滚转无关**；绕枪轴滚转全部由 J6 承担，
而 J6 是 ±360° —— 永远不会成为约束。

### ⛔ 结论：滚转搜索对 J5 无效（此路不通，别做）

"绕枪轴滚转会换解支、从而救 J5" 这个想法在这台机器人上不成立。
θ5 是滚转的不变量。滚转只能改变 θ4/θ6 的分配，而这两个都够宽。
**本方案不含滚转搜索。**（原拟的 C2 已废弃，理由如上。）

### 真正能改 θ5 的量（都在现有代码里，只是作用域不对）

| 量 | 为什么能改 θ5 | 现状 |
|---|---|---|
| **倾角 tilt** | `d` 变化 ⟹ θ5 变化；且 TCP=face+standoff·d 一起移动 | 面组级常量（大平面上恒定） |
| **方位 az** | `d` 在 n 的锥面上绕行 ⟹ θ5 变化 | 行级单值，且白名单卡 ±60° |
| standoff | TCP 沿 d 平移，改变腕心位置 | 工艺量，不轻易动 |
| 基站位姿 | 整体换工作空间 | 需现场摆位 |

`tilt` 和 `az` 都能改 θ5，但**都是全局常量**，这就是大平面上远端被近端绑死的原因。
本方案要解决的就是"作用域"问题。

### 杠杆 C1：倾角里程渐变（bit2，本方案主体）

**驱动量是里程，不是法向**，这是与 bit1 的本质区别。

每行（`computeRows` 定义的一笔）取两个端点：

```
t0 = 起点侧「满足硬约束的最小倾角」   // 网格扫描，不是二分
t1 = 终点侧「满足硬约束的最小倾角」
倾斜量  tilt(s) = t0 + (t1 - t0) * s     s = 行内归一化弧长 ∈ [0,1]
```

- 两端需求相同 ⟹ 退化成面级常量（小面天然向后兼容）
- 两端需求不同 ⟹ 自然产生线性梯度
- 不做"哪个端更差"的判断逻辑，两端各自独立求解即可，简单且不会判错方向

### 杠杆 C1b：方位 az 也按里程渐变（bit2 的可选加强）

同样的毛病：`azRow` 是**行级单值**（`robot_safety.cpp:2010`），一行两端共用一个进刀侧。
而方位对 θ5 的影响与倾角同量级 —— `d` 在 n 的锥面上绕行，θ5 跟着变。

实现上与 C1 完全同构，只是把标量换成 `(tilt, az)` 二元组：

```
端点各求 (t_0, a_0) 与 (t_1, a_1)，中间按里程线性插值
```

**但方位有一处必须先解决**：白名单 `kAzLimitDeg = 60°`（`robot_safety.cpp:2243`）
是工艺边界（"上翻方向左右 60° 内"）。渐变过程中方位可能贴近边界，
需要确认"渐变路径"是否允许擦边，还是要在插值后重新夹一遍。

> 建议：**C1 先只做倾角**。方位渐变引入的工艺争议（枪从侧面斜插）比倾角大得多，
> 等倾角这版跑稳、看过 `rampCappedRows` 的分布再决定要不要扩到方位。

---

## 4. 具体改动（4 个文件）

### 4.1 `robot_safety.h`

```c
#define RS_POSE_FAST        (1 << 0)   // 单次逆解
#define RS_POSE_TILT_PERP   (1 << 1)   // 倾角逐点连续（法向线性）
#define RS_POSE_TILT_RAMP   (1 << 2)   // 【新】倾角沿里程线性渐变
```

`RobotSafetyParams` **尾部**追加 2 个字段（不改动任何既有字段的偏移）：

```c
    // ---- bit2 用：倾角里程渐变 ----
    double rampTiltMaxDeg;   // 渐变可用倾角上限（<=0 时取 maxSprayAngleDeg）
    int    rampProbePct;     // 端点在行内取的位置百分比（默认 8，即 s=0.08 / 0.92）
```

> **ABI 警告（必须写进交接）**：`sizeof(RobotSafetyParams)` 会从 152 变大。
> DLL 与 exe 必须同时重编。既有 `robot_safety_build_stamp()` 已经在界面显示构建戳，
> 用它比对能立刻发现"界面是新的、DLL 是旧的"这种静默失败。
> 若坚持零 ABI 变更，可复用 `upThresholdDeg`（注释已声明"仅保留字段，不参与决策"）
> 承载 `rampTiltMaxDeg` —— 但会让字段语义不可读，**不推荐**。

### 4.2 `robot_safety.cpp`

#### 步骤 1：参数映射

```cpp
p.tiltRamp     = (s.fastPose & RS_POSE_TILT_RAMP) != 0;
p.rampTiltMax  = s.rampTiltMaxDeg > 0 ? s.rampTiltMaxDeg
                                      : p.maxSprayAngleDeg;   // 默认与工艺喷角同口径
p.rampProbePct = (s.rampProbePct >= 0 && s.rampProbePct <= 40) ? s.rampProbePct : 8;
```

`Params` 里加 `bool tiltRamp; double rampTiltMax; int rampProbePct;`

#### 步骤 2：端点探针（新函数，放在 `needTilt` 之后）

沿 `tiltStepDeg` 网格线性扫描，取**第一个满足硬约束**的倾角。
不用二分：`margin` 对倾角不是全局单调的，二分会在非单调段收敛到错的档。

```cpp
// ------------------------------------------------------------------
// 步骤 2：某行端点「满足硬约束的最小倾角」
//
// 为什么按网格线性扫、不二分：jointMargin 沿倾角不是单调函数
// （换解支时会跳），二分依赖单调性，会收敛到错误档位。
// 网格扫的代价可控：0→30° 每 5° 一档 = 7 次 eval，每行只做两次。
// ------------------------------------------------------------------
double probeRowEndTilt(const Evaluator& ev, const Params& p, const TrajIn& t,
                       int i, double az, double cap) {
    V3 nrm = t.normal(i);
    V3 face = t.pos[i] + t.gun[i] * p.standoff;
    V3 axis, uUp;
    tiltAxisAndDir(nrm, p.tiltRefW, axis, uUp);
    if (norm(uUp) < 1e-6) return 0.0;          // 正对基准方向，倾角无意义

    double tMin = 0.0;
    for (double tt = tMin; tt <= cap + 1e-9; tt += p.tiltStepDeg) {
        V3 d = dirOf(nrm, uUp, axis, tt, az);
        PR pr;
        if (!ev.eval(face, d, nullptr, nullptr, tt, az, true, pr)) continue;
        if (!hardOk(pr, p)) continue;
        return tt;                              // 第一个过关的档位
    }
    return -1.0;                                // 封顶仍不满足：这一端做不到
}
```

#### 步骤 3：逐行建剖面（新函数）

```cpp
// ------------------------------------------------------------------
// 步骤 3：按行建立「倾角-里程」线性剖面
//
// 返回 false 表示该行两端都探不出（整行不可达），交给阶段 3 照实报，
// 不在这里自作聪明地给一个"折中倾角"。
// ------------------------------------------------------------------
struct RowRamp {
    bool   ok = false;
    double tilt0 = 0.0, tilt1 = 0.0;   // s=0 / s=1 处的倾角
    int    i0 = -1, i1 = -1;           // 两端实际取的点
};

RowRamp buildRowRamp(const Evaluator& ev, const Params& p, const TrajIn& t,
                     int a, int b, double az) {
    RowRamp r;
    const int len = b - a;
    if (len < 2) return r;
    // 端点往里收一点：行首尾常落在折返圆弧上，那里的法向不代表整行方向
    const int off = std::max(1, (int)(len * p.rampProbePct / 100.0));
    r.i0 = a + off;
    r.i1 = b - off;
    if (r.i1 <= r.i0) { r.i0 = a; r.i1 = b; }

    double t0 = probeRowEndTilt(ev, p, t, r.i0, az, p.rampTiltMax);
    double t1 = probeRowEndTilt(ev, p, t, r.i1, az, p.rampTiltMax);
    if (t0 < 0.0 || t1 < 0.0) return r;        // 至少一端封顶仍不过

    r.ok = true;
    r.tilt0 = t0;
    r.tilt1 = t1;
    return r;
}
```

#### 步骤 4：填 `finalTilt`（修改阶段 2 的分支）

在 `robot_safety.cpp:2028` 那段 `if (p.tiltPerPoint) { ... } else { 面级 }` 里，
**加一个前置分支**，其余一字不动：

```cpp
    if (p.tiltRamp) {
        // ---- 里程渐变口径（RS_POSE_TILT_RAMP）
        // 与面级常量的唯一区别：倾角的作用域从「面组」缩到「行的端点」，中间线性插值。
        // 与逐点连续的区别：bit1 由法向驱动（大平面上恒定），这里由**里程**驱动。
        for (size_t ri = 0; ri < rows.size(); ++ri) {
            const int a = rows[ri].first, b = rows[ri].second;
            const double az = azRow[ri];              // 方位仍取行级决策，不逐点变
            RowRamp rr = buildRowRamp(ev, p, t, a, b, az);
            rampRows[(int)ri] = rr;
            if (!rr.ok) continue;
            const double span = std::max(1, b - a);
            for (int i = a; i <= b; ++i) {
                if (!t.spray[i]) continue;
                const double s = (i - a) / span;      // 归一化里程
                finalTilt[i] = rr.tilt0 + (rr.tilt1 - rr.tilt0) * s;
            }
        }
    } else if (p.tiltPerPoint) {
        ... 原样 ...
    } else {
        ... 原样：面级择优 ...
    }
```

`rampRows` 声明在 `azRow` 旁边，供后面报数用。

#### 步骤 5：阶段 2.5 不要压平梯度

阶段 2.5（`robot_safety.cpp:2571`）的「边界斜坡」会给喷涂段首尾置 `tilt=0`，
再用 `dTiltCap` 双向限速。核对结论：

- 梯度每点增量 = `(t1-t0)/N`。1000mm 行、4mm 点距 = 250 点、差 15° → **0.06°/点**
- `dTiltCap` 在 standoff=100、预算 12mm 时约 **6.9°/点**
- 0.06 ≪ 6.9 ⟹ **梯度不会被限速压掉**，只有首末几点被斜坡归零（与面级常量行为一致）

所以 **阶段 2.5 一行不用改**。这是本方案能低风险落地的主要原因。

> ⚠ 但要把这层关系写进注释：若将来有人把 `kStepBudget` 改到极小，
> 梯度会被限速削平，现象是"开了 ramp 但倾角还是一条水平线"。

#### 步骤 6（可选，先不做）：`eval` 内的滚转搜索 —— 已废弃

**不做。** 理由见第 3 节：θ5 = arccos(d_z)，绕枪轴滚转是 θ5 的不变量，
只能改 θ4/θ6 的分配，而这两个是 ±270°/±360°。对 `JOINT_LIMIT` 零贡献。

留此小节只为记录"已论证过、不要再回头试"，避免下一轮又有人提出同一想法。

### 4.3 `safety_driver.cpp`

```cpp
// -tiltmode 扩成 0|1|2
//   0 = 面级常量（默认）
//   1 = 逐点连续（法向）
//   2 = 【新】里程渐变
if (!strcmp(argv[i], "-tiltmode")) s.fastPose |= atoi(argv[i+1]);
```

另加 `-rampcap <deg>`（默认取 `-sprayang`）。

### 4.4 `wpf/SprayTrajViewer/MainWindow.xaml(.cs)`

现在 `ChkTiltPerPoint`（`MainWindow.xaml:375`）是单个 CheckBox，表达不了三选一。
改成 `ComboBox`「姿态口径」，三项：

```
面级常量（默认） / 逐点连续（法向） / 里程渐变
```

`MainWindow.xaml.cs:2694` 的赋值改为：

```csharp
p.fastPose = (ChkFastPose.IsChecked == true ? 1 : 0)
           | TiltModeBits();      // 由 ComboBox 索引映射：0→0, 1→2, 2→4
```

并同步改 `MainWindow.xaml.cs:2312` 的 `"TiltMode"` 字符串
（现在只有 `per-point` / `face` 两个值，加 `ramp`）。

---

## 5. 报数口径（不做这个，等于改完看不出效果）

`RobotSafetyResult` 已有 `byRule[]`，够用作归因。建议在结论卡加三行：

| 指标 | 含义 | 期望方向 |
|---|---|---|
| `rampMaxStepDeg` | 行内相邻倾角最大步进 | 应是 `(t1-t0)/N` 量级，个位数×0.1 |
| `rampCappedRows` | 封顶仍不满足的行数 | **关键**：>0 说明这是站位/臂展问题 |

**`rampCappedRows > 0` 必须在结论卡直接写明**：

> 该行远端在 30° 喷角上限内仍满足不了关节余量 → `JOINT_LIMIT` 属臂展/站位问题，
> 改倾角无效，需改站位或加导轨。

不写这句，用户会一直以为是算法没调好。

---

## 6. 回归验收（三条硬线）

1. **不开 bit2**：15 个 `cpp/_regress/*.ply` 的 CSV 与 HEAD **逐字节一致**。
   这是项目既有红线，本方案的代码全部走新分支，理论上必然满足，但要实测。
2. **开 bit2 但两端需求相同**（小面/托盘/开口盒）：
   应与默认**逐值一致**（`tilt0==t1` ⟹ 退化为常量）。这是一条很强的自检。
3. **开 bit2 于大平面**：只看四个数 —— 可达率、`byRule[RS_VIOL_JOINT_LIMIT]`、
   `tiltMaxDeg`、行内最大倾角步进。前两个必须降，第三个不得超 30°，第四个必须 <1°。

---

## 7. 风险与边界

| 风险 | 说明 | 缓解 |
|---|---|---|
| 喷角预算是硬上限 | 喷角 ≡ 倾角，远端封顶 30° 仍不过就无解 | 如实报 `rampCappedRows`，不许偷着给 |
| 跨行不连续 | 行与行之间姿态本来就会跳（既有行为），相邻行梯度不同会加大跳幅 | 相邻行 `t0/t1` 差超过 10° 时告警 |
| `TILT_JUMP` 告警可能变多 | 换行处 `t1(行k)` → `t0(行k+1)` 的差可能比面级常量时大 | 与工艺确认换行处可接受跳幅 |
| 探针成本 | 每行 2 端点 × ~7 档 = 14 次 eval | 只在阶段 2 做一次，与 S=400 复核相比可忽略 |
| ABI 变更 | `sizeof(RobotSafetyParams)` 变大 | DLL/exe 必须同编；用 `robot_safety_build_stamp()` 比对 |

---

## 8. 建议实施顺序

1. 只加 bit2 与 `rampTiltMaxDeg` / `rampProbePct` 两个字段，其余一律不动。
2. 先跑回归第 1、2 条（不开位逐字节一致；小面上 `t0==t1` 退化为常量）。
   这两条过了，才说明新分支没污染旧路径。
3. 再跑大平面，只看 `rampCappedRows`：
   - **= 0** ⟹ 分块/渐变有效，继续调 `rampTiltMaxDeg` 找最优
   - **> 0** ⟹ 远端在 30° 内无解，属站位/臂展问题，方案到此为止，
     转去做"站位评估"而不是继续调倾角
4. `rampTiltMaxDeg` 从 15°（半预算）起步，逐步放开到 `maxSprayAngleDeg`。

---

## 9. 需要确认

1. **换行处倾角跳变多大可接受？** 这决定 `rampTiltMaxDeg` 能放多宽。
   相邻行 `t1(行k)` 到 `t0(行k+1)` 的差是既有行为（面级常量时也会跳），
   但渐变会让每行两端不同，跳幅分布会变宽。
2. **是否允许局部超 30° 喷角**换可达性？
   若允许，`rampTiltMaxDeg` 可设 40，余量宽松很多；若不允许，很多大平面会
   直接报 `rampCappedRows > 0`，结论是"需改站位"。
3. 这批次大平面的 **J5 余量具体是负多少**（`marginDeg` 的最小值）。
   如果只差 2~3°，渐变倾角 10° 以内就能救；如果要差 30°，那倾角救不了。

