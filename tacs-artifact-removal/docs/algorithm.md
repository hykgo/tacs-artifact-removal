# 基于先验频率的滑窗最小二乘正弦对消算法：tACS/Hi-tACS 刺激伪迹去除原理

## 摘要

本文介绍一种用于去除 tACS/Hi-tACS 边刺激边采集 EEG/SEEG 中刺激伪迹的算法。该算法基于先验频率假设，将伪迹建模为载频、谐波和边带的窄带正弦叠加，通过多通道平均 PSD 自动识别伪迹频率，构造正弦/余弦参考矩阵，在滑动时间窗内用最小二乘逐段估计伪迹的幅度和相位，并在时域中减除。该方法无需刺激模板，能适应伪迹的缓慢漂移，同时处理多通道信号。本文详细阐述算法的数学原理、实现步骤、参数选择及验证方法。

---
![处理前后 EEG 时域、PSD 与 STFT 对比](../assets/result-comparison.png)
## 1. 引言

tACS（经颅交流电刺激）和 Hi-tACS（高频经颅交流电刺激）在刺激过程中同步采集 EEG/SEEG，可获得刺激对脑活动的调制效应。然而，刺激电流通过电极-皮肤界面在记录信号中引入强伪迹，其幅度通常远大于脑活动，严重遮挡目标信号。

tACS 伪迹的主要特征：

- **窄带正弦**：刺激器输出为固定频率正弦波；
- **多频率成分**：包含载频、谐波（非线性失真）和边带（调制或漂移）；
- **空间一致**：所有通道上频率相同，但幅度和相位不同；
- **非平稳**：幅度和相位随时间缓慢漂移。

本文提出的算法利用这些特征，在已知刺激频率的前提下，通过正弦基最小二乘拟合估计并减除伪迹。

---

## 2. 数学模型

### 2.1 信号模型

将单通道观测信号建模为：

$$
y(t) = s(t) + \sum_{k} \left[ a_k \sin(2\pi f_k t) + b_k \cos(2\pi f_k t) \right] + \varepsilon(t)
$$

其中：

- $y(t)$：观测信号；
- $s(t)$：待保留的脑活动；
- $f_k$：伪迹频率（载频、谐波、边带）；
- $a_k, b_k$：各频率的正弦/余弦系数；
- $\varepsilon(t)$：残余噪声。

### 2.2 正弦的等价表示

任意幅度 $A$、相位 $\phi$ 的正弦可表示为：

$$
A \sin(2\pi f t + \phi) = a \sin(2\pi f t) + b \cos(2\pi f t)
$$

其中：

$$
a = A \cos\phi, \quad b = A \sin\phi
$$

反解：

$$
A = \sqrt{a^2 + b^2}, \quad \phi = \operatorname{atan2}(b, a)
$$

因此，每个频率需要两列基（sin 和 cos）才能覆盖任意幅度和相位。

---

## 3. 算法流程

算法分五个步骤：① 伪迹频率识别 → ② 参考矩阵构造 → ③ 滑窗最小二乘对消 → ④ 背景匹配（可选）→ ⑤ 验证与输出

### 3.1 步骤一：伪迹频率识别

#### 3.1.1 多通道平均 PSD

对所有通道计算 Welch PSD 并求平均：

$$
P_{\text{mean}}(f) = \frac{1}{C} \sum_{c=1}^{C} P_c(f)
$$

多通道平均可增强共同伪迹，抑制随机脑活动。

#### 3.1.2 载频精化

在标称载频 $f_0^{\text{nom}}$ 附近搜索 PSD 最大值，并用抛物线插值精化：

$$
\Delta = \frac{1}{2} \frac{a - g}{a - 2b + g}
$$

真实频率：

$$
f_0 = f_i + \Delta \cdot df
$$

其中 $a, b, g$ 为峰值 bin 及其左右邻居的 PSD 值，$df$ 为频率分辨率。

**抛物线插值原理**：峰值附近近似为抛物线，用三点可确定顶点位置，突破 bin 分辨率限制。

#### 3.1.3 谐波搜索

对 $n = 1, 2, \dots, N_h$，在 $n f_0$ 附近精化并判断显著性：

$$
P_{\text{peak}} > \eta \cdot P_{\text{background}}
$$

其中 $\eta$ 为峰/背景比阈值（默认 10），背景为局部中位数。

#### 3.1.4 边带搜索

对每个谐波 $f_h$，搜索距离在 $[\Delta f_{\min}, \Delta f_{\max}]$ 内且显著的峰，作为边带。

#### 3.1.5 输出

汇总所有频率：

$$
\{f_1, f_2, \dots, f_K\}
$$

### 3.2 步骤二：参考矩阵构造

对每个伪迹频率 $f_k$，构造两列：

$$
\sin(2\pi f_k t), \quad \cos(2\pi f_k t)
$$

将所有列横向拼接：

$$
R = [\sin(2\pi f_1 t), \cos(2\pi f_1 t), \dots, \sin(2\pi f_K t), \cos(2\pi f_K t)]
$$

形状为 $N \times 2K$，$N$ 为采样点数。

**说明**：每个频率两列是必要的，因为任意正弦需要 sin 和 cos 的线性组合才能表示。

### 3.3 步骤三：滑窗最小二乘对消（核心）

#### 3.3.1 全局最小二乘

若假设整段伪迹平稳，可直接解：

$$
\min_{w} \| y - R w \|^2
$$

解析解（正规方程）：

$$
w = (R^T R)^{-1} R^T y
$$

几何意义：将 $y$ 投影到 $R$ 的列空间，投影即伪迹估计。

#### 3.3.2 滑窗的必要性

实际伪迹幅度和相位随时间漂移（电极阻抗变化、刺激器漂移），全局 LS 无法跟踪。因此采用滑窗：每个窗内认为伪迹近似平稳，窗间允许系数变化。

#### 3.3.3 滑窗算法

设定窗长 $W$、步进 $H$（$H < W$）。对每个窗 $[s, s+W)$：

$$
w_j = \arg\min_{w} \| y[s:e] - R[s:e] w \|^2
$$

重建该窗伪迹：

$$
\hat{s}_j[s:e] = R[s:e] w_j
$$

由于窗口重叠，每个采样点被多个窗覆盖。最终伪迹估计为所有覆盖窗的**波形平均**：

$$
\hat{s}(t) = \frac{\sum_j \hat{s}_j(t)}{\sum_j 1}
$$

#### 3.3.4 时域减除

$$
y_{\text{clean}}(t) = y(t) - \hat{s}(t)
$$

**关键点**：平均的是重建波形，不是系数。每个窗的系数只对其覆盖段负责。

#### 3.3.5 数学本质

滑窗 LS 等价于将信号投影到伪迹频率张成的正弦子空间，并逐段剔除该分量。脑活动不在该子空间内，因此被保留在残差中。

### 3.4 步骤四：频域背景匹配（可选）

#### 3.4.1 问题

滑窗 LS 去除伪迹后，PSD 在伪迹频率附近可能出现凹陷（去除了伪迹及其邻近频点）。

#### 3.4.2 方法

对每个伪迹频率 $f_0$：

- 目标区：$|f - f_0| \leq w$；
- 背景区：$w < |f - f_0| \leq w + 3$；
- 背景水平：背景区幅度中位数 $M_{\text{bg}}$。

对目标区每个频点：

$$
|X(f)| < M_{\text{bg}} \Rightarrow |X(f)| \leftarrow M_{\text{bg}}
$$

相位保持不变。

#### 3.4.3 注意事项

- 此步骤人为补功率，仅用于改善 PSD 外观；
- **做功率分析、连通性分析时必须关闭**，否则引入虚假功率；
- 相位不变可避免时域波形突变。

### 3.5 步骤五：验证与输出

**残差报告**：比较处理前后 PSD，计算伪迹频率处功率下降百分比。

**可视化**：

- 时域波形对比：确认伪迹幅度降低；
- PSD 对比：确认目标频率下降、非目标频段不变；
- STFT 时频图：确认伪迹频率处能量被有效抑制。

**输出**：保存为 FIF/EDF 格式。

---

## 4. 参数选择

| 参数 | 含义 | 建议 |
|---|---|---|
| $f_0^{\text{nom}}$ | 标称载频 | 必须准确 |
| $N_h$ | 谐波数 | 8 左右 |
| $\eta$ | 峰/背景比阈值 | 10 左右 |
| $W$ | 窗长 | 1 s |
| $H$ | 步进 | 0.1 s |
| $w$ | 背景匹配宽度 | 1.5 Hz |
| 背景匹配 | 是否补功率 | 功率分析时关闭 |

---

## 5. 算法特性

### 5.1 优点

- 无需刺激模板或触发标记；
- 同时处理载频、谐波、边带；
- 每通道独立，适应不同电极的幅相差异；
- 滑窗适应伪迹缓慢漂移；
- 数学原理清晰，可解释性强。

### 5.2 局限

- 假设伪迹为稳定窄带正弦；
- 与脑活动同频时会误删；
- 背景匹配是双刃剑；
- 频率搜索依赖参数；
- 不处理刺激开/关瞬态。

### 5.3 注意事项

- 边带搜索范围应避开 delta 频段；
- 频率过多或窗过短会导致过拟合；
- 刺激瞬态需分段处理或额外去漂移；
- 建议保留原始数据以便回溯。

---

## 6. 总结

本文介绍的算法通过“PSD 找频率 → sin/cos 建参考 → 滑窗最小二乘估幅相 → 时域减除 → 验证输出”的流程，有效去除 tACS/Hi-tACS 刺激伪迹。其核心是将信号投影到伪迹频率张成的正弦子空间并剔除该分量，滑窗机制适应伪迹的非平稳变化。该算法原理直观、实现简洁，适用于已知刺激频率的 EEG/SEEG 伪迹去除。

---

**关键词**：tACS；Hi-tACS；刺激伪迹去除；最小二乘；滑窗；正弦对消；EEG；SEEG



```
import numpy as np
import matplotlib.pyplot as plt
import mne
from scipy.signal import welch, find_peaks, stft

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ============================================================
# 1. 读 EDF
# ============================================================
file = r"R2_Val40_During.edf"

raw = mne.io.read_raw_edf(file, preload=True, verbose='error')
fs = raw.info['sfreq']
print(f"采样率: {fs} Hz, 通道数: {len(raw.ch_names)}, 时长: {raw.times[-1]:.2f} s")

data = raw.get_data()
n_ch, n_samp = data.shape
t = raw.times

# ============================================================
# 2. 多通道平均 PSD
# ============================================================
f_psd, Pxx_all = welch(data, fs=fs,
                       nperseg=int(fs*8), noverlap=int(fs*4), axis=-1)
Pxx_mean = Pxx_all.mean(axis=0)

# ============================================================
# 3. 精化频率 + 自动找边带
# ============================================================
def refine_peak_near(f_psd, P, f_target, search=2.0):
    mask = (f_psd >= f_target - search) & (f_psd <= f_target + search)
    if mask.sum() == 0:
        return f_target
    idx_local = np.argmax(P[mask])
    i = np.where(mask)[0][idx_local]
    if i <= 0 or i >= len(P) - 1:
        return f_psd[i]
    a, b, g = P[i-1], P[i], P[i+1]
    denom = a - 2*b + g
    if denom == 0:
        return f_psd[i]
    d = 0.5 * (a - g) / denom
    return f_psd[i] + d * (f_psd[1] - f_psd[0])

f0_nominal = 40.0
f0_real = refine_peak_near(f_psd, Pxx_mean, f0_nominal, search=2.0)
print(f"载频: 标称 {f0_nominal} → 实测 {f0_real:.4f} Hz")

peaks, _ = find_peaks(Pxx_mean,
                      height=np.median(Pxx_mean) * 10,
                      distance=max(1, int(1.0 / (f_psd[1]-f_psd[0]))))

harmonics = []
for n in range(1, 8):
    f_t = n * f0_real
    if f_t >= fs / 2:
        break
    f_real = refine_peak_near(f_psd, Pxx_mean, f_t, search=2.0)
    p_peak = Pxx_mean[np.argmin(np.abs(f_psd - f_real))]
    p_bg = np.median(Pxx_mean[(f_psd > f_real - 5) & (f_psd < f_real + 5)])
    if p_peak > p_bg * 10:
        harmonics.append(f_real)
        print(f"  谐波 {f_t:6.2f} → {f_real:8.4f} Hz  (峰/背景 = {p_peak/p_bg:.1f})")

sidebands = []
for f_h in harmonics:
    for f_pk in f_psd[peaks]:
        if f_h == f_pk:
            continue
        if 2.0 < abs(f_pk - f_h) < 30.0:
            p_peak = Pxx_mean[np.argmin(np.abs(f_psd - f_pk))]
            p_bg = np.median(Pxx_mean[(f_psd > f_pk - 3) & (f_psd < f_pk + 3)])
            if p_peak > p_bg * 10:
                sidebands.append(f_pk)

sidebands = sorted(set(round(f, 4) for f in sidebands))
if sidebands:
    print(f"检测到边带: {sidebands}")

all_freqs = sorted(set([round(f, 4) for f in harmonics] + sidebands))
print(f"\n最终使用伪迹频率: {all_freqs}")

plt.figure(figsize=(10, 4))
plt.semilogy(f_psd, Pxx_mean + 1e-30)
for f in all_freqs:
    plt.axvline(f, color='r', ls='--', alpha=0.5)
plt.xlim(0, min(200, fs/2))
plt.xlabel("Hz"); plt.ylabel("PSD")
plt.title("多通道平均 PSD 与选中的伪迹频率")
plt.tight_layout(); plt.show()

# ============================================================
# 4. 滑窗 LS 对消
# ============================================================
def build_ref(t, freqs):
    cols = []
    for f in freqs:
        cols.append(np.sin(2 * np.pi * f * t))
        cols.append(np.cos(2 * np.pi * f * t))
    return np.column_stack(cols)

R = build_ref(t, all_freqs)
print(f"\n参考矩阵形状: {R.shape}")

def sliding_ls_cancel(y, R, win, hop):
    n = len(y)
    est_sum = np.zeros(n)
    count   = np.zeros(n)

    for s in range(0, n - win + 1, hop):
        e = s + win
        ws, *_ = np.linalg.lstsq(R[s:e], y[s:e], rcond=None)
        est_sum[s:e] += R[s:e] @ ws
        count[s:e]   += 1

    if (n - win) % hop != 0:
        s = n - win
        ws, *_ = np.linalg.lstsq(R[s:], y[s:], rcond=None)
        est_sum[s:] += R[s:] @ ws
        count[s:]   += 1

    count[count == 0] = 1
    return y - est_sum / count

win = int(fs * 1.0)
hop = int(fs * 0.1)

print(f"开始滑窗 LS（win={win/fs:.1f}s, hop={hop/fs:.2f}s）...")
data_clean = np.zeros_like(data)
for ch in range(n_ch):
    data_clean[ch] = sliding_ls_cancel(data[ch], R, win, hop)
    if (ch + 1) % 10 == 0:
        print(f"  已处理 {ch+1}/{n_ch} 通道")

# ============================================================
# 5. 频域背景匹配（覆盖整个凹陷宽度）
# ============================================================
print("\n开始频域背景匹配...")

def background_match_fft(x, fs, freqs, half_width=1.5):
    """
    把 freqs 附近 ±half_width Hz 内、幅度低于邻近背景中位数的所有 bin
    抬到背景中位数，相位不变。
    """
    X = np.fft.rfft(x)
    freqs_fft = np.fft.rfftfreq(len(x), 1/fs)

    for f0 in freqs:
        # 需要修补的目标区域
        mask_target = np.abs(freqs_fft - f0) <= half_width
        # 用于估计背景的 bin：目标 ±half_width 之外，±(half_width+3) 之内
        mask_bg = (np.abs(freqs_fft - f0) > half_width) & \
                  (np.abs(freqs_fft - f0) <= half_width + 3.0)
        if mask_bg.sum() < 10:
            continue

        bg_median = np.median(np.abs(X[mask_bg]))

        idx_target = np.where(mask_target)[0]
        for i in idx_target:
            if np.abs(X[i]) < bg_median:
                X[i] = bg_median * np.exp(1j * np.angle(X[i]))

    return np.fft.irfft(X, n=len(x))

for ch in range(n_ch):
    data_clean[ch] = background_match_fft(data_clean[ch], fs, all_freqs,
                                          half_width=1.5)
    if (ch + 1) % 10 == 0:
        print(f"  背景匹配 {ch+1}/{n_ch} 通道")

# ============================================================
# 6. 检查残差
# ============================================================
f2, P2_all = welch(data_clean, fs=fs,
                   nperseg=int(fs*8), noverlap=int(fs*4), axis=-1)
P2_mean = P2_all.mean(axis=0)

print("\n处理后残差（多通道平均）：")
for f in all_freqs:
    p_before = Pxx_mean[np.argmin(np.abs(f_psd - f))]
    p_after  = P2_mean[np.argmin(np.abs(f_psd - f))]
    print(f"  {f:8.4f} Hz : {p_before:.3e} → {p_after:.3e}  "
          f"(降到 {p_after/p_before*100:.2f}%)")

# ============================================================
# 7. 可视化
# ============================================================
ch_show = 12
print(f"\n展示通道: {raw.ch_names[ch_show]}")

seg = slice(int(2*fs), int(2.5*fs))
fig, axes = plt.subplots(3, 2, figsize=(14, 10))

axes[0, 0].plot(t[seg], data[ch_show, seg], label="原始", lw=0.8)
axes[0, 0].plot(t[seg], data_clean[ch_show, seg], label="处理后", lw=0.8)
axes[0, 0].set_title(f"时域波形 ({raw.ch_names[ch_show]}, 2.0-2.5 s)")
axes[0, 0].legend(); axes[0, 0].set_xlabel("t (s)")

axes[0, 1].plot(t, data[ch_show], lw=0.5, alpha=0.7, label="原始")
axes[0, 1].plot(t, data_clean[ch_show], lw=0.5, alpha=0.7, label="处理后")
axes[0, 1].set_title("全段波形")
axes[0, 1].legend()

f1, P1 = welch(data[ch_show],        fs=fs, nperseg=int(fs*8), noverlap=int(fs*4))
f2c, P2c = welch(data_clean[ch_show], fs=fs, nperseg=int(fs*8), noverlap=int(fs*4))

axes[1, 0].semilogy(f1, P1 + 1e-30, label="原始", lw=1)
axes[1, 0].semilogy(f2c, P2c + 1e-30, label="处理后", lw=1)
for f in all_freqs:
    axes[1, 0].axvline(f, color='r', ls='--', alpha=0.3)
axes[1, 0].set_xlim(0, min(200, fs/2))
axes[1, 0].set_title("PSD 对比 (dB)")
axes[1, 0].legend(); axes[1, 0].set_xlabel("Hz"); axes[1, 0].set_ylabel("PSD")

f_lo = max(0, f0_real - 20); f_hi = min(fs/2, f0_real + 20)
axes[1, 1].semilogy(f1, P1 + 1e-30, label="原始", lw=1)
axes[1, 1].semilogy(f2c, P2c + 1e-30, label="处理后", lw=1)
axes[1, 1].set_xlim(f_lo, f_hi)
axes[1, 1].set_title(f"PSD 放大 ({f_lo:.0f}-{f_hi:.0f} Hz, dB)")
axes[1, 1].legend(); axes[1, 1].set_ylabel("PSD")

def plot_stft(ax, x, fs, title):
    f, tt, Z = stft(x, fs=fs, nperseg=int(fs*0.5), noverlap=int(fs*0.4))
    ax.pcolormesh(tt, f, 20*np.log10(np.abs(Z)+1e-12),
                  shading='gouraud', cmap='magma')
    ax.set_ylim(0, min(200, fs/2))
    ax.set_title(title); ax.set_xlabel("t (s)"); ax.set_ylabel("Hz")

plot_stft(axes[2, 0], data[ch_show],       fs, "STFT 原始 (dB)")
plot_stft(axes[2, 1], data_clean[ch_show], fs, "STFT 处理后 (dB)")

plt.tight_layout(); plt.show()

# ============================================================
# 8. 保存
# ============================================================
out_edf = file.replace(".edf", "_clean.edf")
out_fif = file.replace(".edf", "_clean_raw.fif")

raw_clean = mne.io.RawArray(data_clean, raw.info.copy(), verbose='error')
raw_clean.save(out_fif, overwrite=True)
print(f"\n已保存 FIF: {out_fif}")

try:
    mne.export.export_raw(out_edf, raw_clean, fmt='edf', overwrite=True)
    print(f"已保存 EDF: {out_edf}")
except Exception as e:
    print(f"EDF 保存失败: {e}")
```
