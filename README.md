# CERTEO: scene-level certification of SAR-to-optical translation

This repository contains the code for the paper *Scene-Level Certification of SAR-to-Optical Translation: The Certification Tax and Capped Trust Masks* (under review).

A SAR-to-optical generator turns a radar patch into an image that looks optical, but it does not say which of its pixels can be trusted. CERTEO builds a **trust mask** for any frozen generator, with a stated guarantee on the share of pixels that are trusted but wrong. This README explains every piece of mathematics the code uses and how to run the code. The same material, typeset as a paper-style document, is in [`docs/certeo_method.pdf`](docs/certeo_method.pdf); its LaTeX source is in [`docs/certeo_method.tex`](docs/certeo_method.tex).

**Main findings**
1. Calibrating on patches, as Earth-observation benchmarks usually do, is *exactly* scene-level calibration at a higher level than the one requested. In our experiments it left 41%–51% of deployed users over their error budget, against 0%–20% for scene-level calibration.
2. Every distribution-free per-user certificate has to pay a **certification tax** proportional to $`(q-\alpha)/m`$, where $`m`$ is the number of labelled scenes and $`q`$ bounds the loss. With $`q = 1`$, no method can certify anything from fewer than 45 scenes.
3. A **capped trust mask** (trust at most a share $`q`$ of every patch) lowers the tax. Combined with a closed-form **product certificate**, it certifies 7%–11% of pixels from 8–16 scenes, where the uncapped exact bound certifies none.

---

## Contents

1. [Setting and notation](#1-setting-and-notation)
2. [Conformal risk control with scenes as units](#2-conformal-risk-control-with-scenes-as-units)
3. [Proposition 1: unit equivalence](#3-proposition-1-unit-equivalence)
4. [Proposition 2: how often a user is over budget](#4-proposition-2-how-often-a-user-is-over-budget)
5. [Proposition 3: the certification tax](#5-proposition-3-the-certification-tax)
6. [Capped trust masks](#6-capped-trust-masks)
7. [Proposition 4: the product certificate](#7-proposition-4-the-product-certificate)
8. [The cap rule](#8-the-cap-rule)
9. [The full procedure](#9-the-full-procedure)
10. [Identifiability floor and ceiling](#10-identifiability-floor-and-ceiling)
11. [A certificate on the error among trusted pixels](#11-a-certificate-on-the-error-among-trusted-pixels)
12. [Computing everything from histograms](#12-computing-everything-from-histograms)
13. [Experimental design](#13-experimental-design)
14. [Results](#14-results)
15. [What is certified, and what is not](#15-what-is-certified-and-what-is-not)
16. [Running the notebook](#16-running-the-notebook)
17. [Repository layout](#17-repository-layout)
18. [References](#18-references)

---

## 1. Setting and notation

| Symbol | Meaning |
|---|---|
| $`x`$, $`y`$, $`\hat y`$ | SAR patch, reference optical patch, generated optical patch |
| $`b \in \{1,2,3\}`$ | optical band (red, green, blue), scaled to $`[0,1]`$ |
| $`\tau`$ | error tolerance: a pixel is **wrong** if $`\max_b \lvert \hat y_b - y_b \rvert > \tau`$ |
| $`s`$ | trust score per pixel; higher means *less* trustworthy |
| $`\lambda`$ | threshold; the trust mask is $`\{ s \le \lambda \}`$ |
| $`P`$ | number of pixels in a patch ($`256 \times 256`$) |
| $`n`$ | patches per scene |
| $`m`$ | number of labelled calibration scenes |
| $`\alpha`$ | error budget (target false-trust share) |
| $`\delta`$ | allowed probability of a failed certificate |
| $`q`$ | cap: the largest share of a patch that may be trusted |

A **scene** is one region in one season. Patches from one scene share the weather, season, land cover and sensor state, so they are strongly correlated.

The **loss of a patch** is the share of its pixels that are trusted *and* wrong:

```math
\ell(\lambda) = \frac{1}{P} \sum_{p=1}^{P} \mathbf{1}\{ s_p \le \lambda \} \, \mathbf{1}\Big\{ \max_b \lvert \hat y_{p,b} - y_{p,b} \rvert > \tau \Big\}.
```

The **loss of scene** $`i`$ is the mean over its $`n`$ patches, $`L_i(\lambda) = \frac{1}{n}\sum_{k=1}^{n} \ell_{ik}(\lambda)`$. It never decreases as $`\lambda`$ grows, because a larger threshold only adds trusted pixels.

The **risk** of a threshold is the expected scene loss, $`R(\lambda) = \mathbb{E}[L_i(\lambda)]`$. The **certified area** is the share of pixels that are trusted at the selected threshold.

**Worked example.** A patch has $`P = 100`$ pixels. With threshold $`\lambda`$, 30 pixels are trusted, and 4 of these 30 are wrong. Then $`\ell(\lambda) = 4/100 = 0.04`$: the loss counts wrong pixels over *all* pixels, not over the trusted ones. The error among trusted pixels would be $`4/30 \approx 0.13`$; that is a different quantity, discussed in Section 15.

---

## 2. Conformal risk control with scenes as units

Conformal risk control (CRC) [1] chooses a threshold from labelled calibration units. Take $`m`$ exchangeable scenes whose losses lie in $`[0, B]`$. Scene-level CRC selects

```math
\hat\lambda = \max \Big\{ \lambda : \frac{\sum_{i=1}^{m} L_i(\lambda) + B}{m+1} \le \alpha \Big\}.
```

The guarantee is on a new scene, on average over calibration sets:

```math
\mathbb{E}\big[ L_{m+1}(\hat\lambda) \big] \le \alpha .
```

The "$`+B`$" in the numerator is the price of not knowing the new scene's loss. It also means CRC needs a minimum number of scenes: even with zero observed loss, it accepts only if $`B/(m+1) \le \alpha`$, that is, $`m \ge B/\alpha - 1`$. With $`B = 1`$ and $`\alpha = 0.05`$ this is 19 scenes.

We say a user is **over budget** if, for the calibration set this user happened to draw, $`R(\hat\lambda) > \alpha`$. CRC controls the *average* risk, not this per-user event.

---

## 3. Proposition 1: unit equivalence

**Patch-level CRC** applies the rule of Section 2 to the $`mn`$ patches as if they were exchangeable, with $`B = 1`$.

> **Proposition 1.** With $`n`$ patches in each of $`m`$ scenes and $`B = 1`$, patch-level CRC at level $`\alpha`$ selects exactly the same threshold as scene-level CRC at level
>
> ```math
> \alpha' = \alpha + \frac{(1-\alpha)(1 - 1/n)}{m+1}.
> ```

**Proof.** The sum over patches is $`\sum_{i,k} \ell_{ik} = n \sum_i L_i`$. Patch-level CRC accepts $`\lambda`$ if

```math
\frac{n \sum_i L_i(\lambda) + 1}{mn + 1} \le \alpha
\iff
\sum_i L_i(\lambda) \le \frac{\alpha (mn+1) - 1}{n}.
```

Scene-level CRC at level $`\alpha'`$ accepts $`\lambda`$ if $`\sum_i L_i(\lambda) \le \alpha'(m+1) - 1`$. Setting the two right-hand sides equal gives

```math
\alpha' = \frac{\alpha(mn+1) - 1 + n}{n(m+1)} = \alpha + \frac{(1-\alpha)(1-1/n)}{m+1}.
```

The two rules accept the same thresholds, so they select the same $`\hat\lambda`$. Since scenes, not patches, are exchangeable, the CRC theorem gives only $`\mathbb{E}[L_{m+1}(\hat\lambda)] \le \alpha'`$. $`\square`$

**What it means.** This sharpens the observation of Dunn et al. [10] that pooling grouped data gives only approximate guarantees. Cutting more patches from each scene does not help: as $`n`$ grows, $`\alpha'`$ approaches $`\alpha + (1-\alpha)/(m+1)`$, not $`\alpha`$. Example: $`m = 32`$, $`n = 64`$, $`\alpha = 0.05`$ gives $`\alpha' = 0.078`$.

**Check on real data.** In 216 configurations of generator, score, budget, $`m \in \{16, 32, 64\}`$ and $`n \in \{1, 4, 16, 64\}`$ (400 draws each), the two selections never differed. The patch-level mean risk exceeded $`\alpha`$ in 47 configurations, and never exceeded $`\alpha'`$.

---

## 4. Proposition 2: how often a user is over budget

Let $`\lambda^{+}`$ be the first threshold on the grid with $`R(\lambda^{+}) > \alpha`$, and let $`\bar L = \frac{1}{m}\sum_i L_i(\lambda^{+})`$ be its calibration mean. Because losses are monotone in $`\lambda`$, the user is over budget **exactly** when CRC accepts $`\lambda^{+}`$, which happens exactly when

```math
\bar L \le t, \qquad
t = \begin{cases}
\alpha - \dfrac{1-\alpha}{m} & \text{scene-level CRC},\\[2ex]
\alpha - \dfrac{1-\alpha}{mn} & \text{patch-level CRC}.
\end{cases}
```

> **Proposition 2 (skew-corrected law).** Let $`\sigma_n^2`$ and $`\gamma`$ be the variance and skewness of $`L_i(\lambda^{+})`$, and $`z = (t - R(\lambda^{+}))\sqrt{m}/\sigma_n`$. Then
>
> ```math
> \mathbb{P}(\text{over budget}) \approx \Phi(z) - \frac{\gamma \, (z^2 - 1) \, \phi(z)}{6 \sqrt{m}},
> ```
>
> where $`\Phi`$ and $`\phi`$ are the standard normal distribution function and density.

This is the one-term Edgeworth expansion [2] of the distribution of the mean $`\bar L`$.

**Why patch-level calibration fails about half of the users.** The variance of a scene loss splits into a between-scene and a within-scene part:

```math
\sigma_n^2 = \sigma_b^2 + \frac{\sigma_w^2}{n}.
```

The between-scene part $`\sigma_b^2`$ stays however many patches are cut. Patch-level CRC sets $`t`$ almost equal to $`\alpha`$, so $`t - R(\lambda^{+})`$ shrinks to zero, and $`z \to 0`$. Then

```math
\mathbb{P}(\text{over budget}) \to \frac{1}{2} + \frac{\gamma}{6\sqrt{2\pi m}},
```

which is *above* one half when scene losses are right-skewed ($`\gamma > 0`$), as they are on real scenes (skewness 2.5–4.2).

**Check on real data.** Over 126 configurations (three generators, $`m \in \{16, 32, 64\}`$, $`n \in \{1, \dots, 64\}`$, 2,000 draws), the largest error of the law is 0.031 (mean 0.006), against 0.082 (mean 0.021) for the plain normal approximation $`\Phi(z)`$.

---

## 5. Proposition 3: the certification tax

A **per-user certificate** asks for more than CRC:

```math
\mathbb{P}\big( R(\hat\lambda) \le \alpha \big) \ge 1 - \delta ,
```

where the probability is over the calibration set. A procedure is $`(\alpha, \delta)`$-valid if this holds for every distribution of scene losses on $`[0, q]`$.

> **Proposition 3 (certification tax).** Let scene losses be i.i.d. in $`[0, q]`$. Any $`(\alpha,\delta)`$-valid procedure accepts a threshold of risk $`R < \alpha`$ with probability at most
>
> ```math
> \delta \left( \frac{q - R}{q - \alpha} \right)^{m}.
> ```
>
> Accepting it with probability at least $`\beta`$ therefore requires
>
> ```math
> R \le \alpha - (q - \alpha)\Big[ (\beta/\delta)^{1/m} - 1 \Big].
> ```

**Proof** (in the spirit of Bahadur and Savage [9]). Let $`P_0`$ be a loss distribution on $`[0,q]`$ with mean $`R < \alpha`$, and $`A`$ the event that the procedure accepts. Mix in a point mass at $`q`$:

```math
P_\varepsilon = (1-\varepsilon) P_0 + \varepsilon \, \delta_q , \qquad \varepsilon \in (0,1).
```

Its mean is $`(1-\varepsilon)R + \varepsilon q`$, which exceeds $`\alpha`$ once $`\varepsilon > (\alpha - R)/(q - R)`$, so validity demands $`P_\varepsilon^m(A) \le \delta`$. With probability $`(1-\varepsilon)^m`$ all $`m`$ scenes come from the $`P_0`$ part, so $`P_\varepsilon^m(A) \ge (1-\varepsilon)^m P_0^m(A)`$ and

```math
P_0^m(A) \le \frac{\delta}{(1-\varepsilon)^m} \;\longrightarrow\; \delta \left( \frac{q-R}{q-\alpha} \right)^m
\quad \text{as } \varepsilon \downarrow \frac{\alpha - R}{q - R}.
```

Requiring $`P_0^m(A) \ge \beta`$ and solving for $`R`$ gives the inequality. $`\square`$

The subtracted term $`(q-\alpha)[(\beta/\delta)^{1/m} - 1]`$ is the **tax**; for large $`m`$ it is about $`(q-\alpha)\ln(\beta/\delta)/m`$. CRC pays the similar term $`(q-\alpha)/m`$. Both are proportional to $`(q - \alpha)/m`$.

**Minimum number of scenes.** With $`R = 0`$ and $`\beta = 1`$, nothing can be accepted unless

```math
m \ge \frac{\ln(1/\delta)}{\ln\big( q / (q - \alpha) \big)} .
```

At $`\alpha = 0.05`$, $`\delta = 0.1`$:

| Cap $`q`$ | 1 | 0.5 | 0.25 | 0.15 |
|---|---|---|---|---|
| CRC (average guarantee), $`m \ge q/\alpha - 1`$ | 19 | 9 | 4 | 2 |
| Any per-user certificate (Proposition 3) | 45 | 22 | 11 | 6 |
| Product certificate, $`c = 0.75`$ (Section 7) | 60 | 29 | 14 | 8 |

**Intuition for 45.** If 5.1% of scenes were completely wrong, a sample of 44 scenes would contain none of them with probability $`0.949^{44} \approx 0.10`$. Such a sample looks exactly like a perfect world, so no honest method may certify from it.

**Check.** In simulation, the product certificate with $`c = 1`$ accepts an all-zero sample from exactly these values of $`m`$, and the mixture adversary is then accepted with probability 0.089–0.098, so the bound is attained.

---

## 6. Capped trust masks

The cap $`q`$ is the only term of the tax that the user controls. A **capped trust mask** trusts the pixels with $`s \le \lambda`$, but at most a share $`q`$ of each patch, taking the lowest scores first:

```math
M_q(\lambda) = \{ s \le \lambda \} \cap \{ s \le s_{(\lceil qP \rceil)} \},
```

where $`s_{(k)}`$ is the $`k`$-th smallest score in the patch. Every patch loss, and hence every scene loss, then lies in $`[0, q]`$, and the loss is still nondecreasing in $`\lambda`$.

**Example.** At $`\alpha = 0.05`$, $`\delta = 0.1`$, moving from $`q = 1`$ to $`q = 0.25`$ lowers the minimum number of scenes for any per-user certificate from 45 to 11.

**The trade-off.** The certified area can never exceed $`q`$. A small cap costs area; a large cap costs tax. Section 8 balances the two.

---

## 7. Proposition 4: the product certificate

> **Proposition 4.** For i.i.d. scene losses in $`[0, q]`$, $`q > \alpha`$, and a fixed bet $`c \in (0, 1]`$, let $`\hat\lambda`$ be the largest $`\lambda`$ with
>
> ```math
> \sum_{i=1}^{m} \log\Big( 1 + c \, \frac{\alpha - L_i(\lambda)}{q - \alpha} \Big) \ge \log \frac{1}{\delta} .
> ```
>
> Then $`\mathbb{P}\big( R(\hat\lambda) > \alpha \big) \le \delta`$.

**Proof.** At the first unsafe threshold $`\lambda^{+}`$, let

```math
W_i = 1 + c \, \frac{\alpha - L_i(\lambda^{+})}{q - \alpha}.
```

Since $`L_i \le q`$ and $`c \le 1`$, $`W_i \ge 1 - c \ge 0`$. Its mean is

```math
\mathbb{E}[W_i] = 1 + c \, \frac{\alpha - R(\lambda^{+})}{q - \alpha} < 1 .
```

By independence $`\mathbb{E}\big[\prod_i W_i\big] \le 1`$, and Markov's inequality gives $`\mathbb{P}\big(\prod_i W_i \ge 1/\delta\big) \le \delta`$. The statistic never increases with $`\lambda`$, because every $`L_i`$ never decreases, so the event $`\{\hat\lambda \ge \lambda^{+}\}`$ is exactly the event $`\{\prod_i W_i \ge 1/\delta\}`$. Risk is nondecreasing in $`\lambda`$, so $`R(\hat\lambda) > \alpha`$ implies $`\hat\lambda \ge \lambda^{+}`$. $`\square`$

**Interpretation as a bet.** Each scene is a bet that pays $`1 + c(\alpha - L_i)/(q - \alpha)`$ on a stake of 1. If the threshold were unsafe, the bet would lose money on average, so the chance of multiplying the stake to $`1/\delta`$ is at most $`\delta`$. The certificate is a fixed-bet member of the betting family of Waudby-Smith and Ramdas [3]. It is closed-form, and it does not depend on the order of the scenes.

**Worked example.** $`m = 16`$, $`q = 0.13`$, $`\alpha = 0.05`$, $`c = 0.75`$, $`\delta = 0.1`$. A scene with loss 0.02 adds

```math
\log\Big(1 + 0.75 \cdot \frac{0.05 - 0.02}{0.13 - 0.05}\Big) = \log(1.28125) \approx 0.248 .
```

Sixteen such scenes give about 3.96, which exceeds $`\log(1/0.1) \approx 2.30`$, so the threshold is certified.

**Relation to Proposition 3.** If all observed losses equal $`r`$, the certificate accepts when

```math
r \le \alpha - \frac{(q - \alpha)\big(\delta^{-1/m} - 1\big)}{c},
```

which matches the necessary condition of Proposition 3 up to the factor $`1/c`$ (and $`\beta = 1`$).

---

## 8. The cap rule

We spend a quarter of the budget on the tax. Let $`u_m = \delta^{-1/m} - 1`$ for the product certificate and $`u_m = 1/m`$ for CRC. The cap is

```math
q^{\star}(m) = \min\Big\{ 1, \; \alpha \big[ 1 + 0.25 / u_m \big] \Big\}.
```

**Where it comes from.** The tax is about $`(q - \alpha)\,u_m`$. Setting it to $`0.25\,\alpha`$ and solving for $`q`$ gives $`q = \alpha(1 + 0.25/u_m)`$. As $`m`$ grows, $`u_m \to 0`$, so $`q^\star(m) \to 1`$: with many scenes, no cap is needed.

**Example.** $`m = 32`$, $`\delta = 0.1`$, $`\alpha = 0.05`$: $`u_{32} = 10^{1/32} - 1 \approx 0.0746`$, so $`q^{\star} \approx 0.05 \times (1 + 3.35) \approx 0.218`$.

We fix $`c = 0.75`$. The constants 0.25 and 0.75 were chosen on the calibration scenes; the held-out winter scenes then tested them without further tuning.

---

## 9. The full procedure

**Algorithm 1: scene-level certification with a capped trust mask.**

Input: frozen generator, score $`s`$, $`m`$ labelled scenes, budget $`\alpha`$, confidence $`\delta`$, tolerance $`\tau`$.

1. Set the cap $`q \leftarrow q^{\star}(m)`$.
2. For each scene $`i`$ and each threshold $`\lambda`$ on a fixed grid:
   1. trust the pixels with $`s \le \lambda`$, at most a share $`q`$ of each patch;
   2. compute $`L_i(\lambda)`$, the mean share of trusted pixels with error above $`\tau`$.
3. Select $`\hat\lambda`$:
   - per-user guarantee: the largest $`\lambda`$ that satisfies the product certificate (Section 7);
   - average guarantee: CRC (Section 2) with $`B = q`$.
4. If no threshold qualifies, report the minimum $`m`$ from Proposition 3 and refuse.

Output: the capped mask at $`\hat\lambda`$, applied to every new scene.

---

## 10. Identifiability floor and ceiling

Some error is not the generator's fault: where the radar does not determine the optical image, *no* generator can be right reliably. For a generator that sees only the SAR patch, the output is independent of the true optical value given the radar. Hence, for every pixel,

```math
\mathbb{P}\big( \lvert \hat y_b - y_b \rvert > \tau \mid x \big) \ge 1 - Q_{\tau,b}(x),
\qquad
Q_{\tau,b}(x) = \sup_{c} \, \mathbb{P}\big( \lvert y_b - c \rvert \le \tau \mid x \big).
```

$`Q_{\tau,b}(x)`$ is the largest conditional probability that any window of width $`2\tau`$ can hold. The **floor** of a pixel is $`f = \max_b \big[1 - Q_{\tau,b}(x)\big]`$.

**Example.** If $`p(y \mid x)`$ puts half of its mass at 0.2 and half at 0.8, and $`\tau = 0.1`$, any guess is within $`\tau`$ of at most one of the two values, so the error rate at that pixel is at least 50%.

**Ceiling.** Every trust rule pays at least the floor for each pixel it trusts. So the largest area any rule can certify at budget $`\alpha`$ comes from filling pixels in order of increasing floor until the summed floor reaches $`\alpha`$:

```math
\mathrm{CTA}^{\star}(\alpha) = \max_{S} \frac{\mathbb{E}\lvert S \rvert}{P}
\quad \text{subject to} \quad
\mathbb{E}\sum_{p \in S} f_p \le \alpha P .
```

This is a Chow-type rejection bound [4]. We estimate $`f`$ with a per-pixel mixture density network (four Gaussian components per band).

**Caveat.** The estimated floor is a *valid but loose* lower bound: in the lowest floor decile it predicts a miss rate of 0.034, while the $`\tau`$-mode generator built from the same network misses 0.289. The resulting ceiling (0.441 at $`\alpha = 0.05`$) is therefore an overestimate.

---

## 11. A certificate on the error among trusted pixels

Users may prefer "of what I trust, at most $`\epsilon`$ is wrong". With $`K(\lambda)`$ the trusted share, this is the hypothesis

```math
\mathbb{E}[L(\lambda)] \le \epsilon \, \mathbb{E}[K(\lambda)] .
```

We test it with the product test applied to $`L - \epsilon K \in [-\epsilon, 1 - \epsilon]`$, with a Bonferroni correction over nine fixed thresholds (the ratio is not monotone in $`\lambda`$). With $`\epsilon = 0.3`$, $`\tau = 0.2`$ and 64 scenes it certifies 63% (Pix2Pix), 75% (U-Net) and 100% (cloud removal) of pixels; with 32 scenes, 8%, 11% and 84%. At $`\epsilon = 0.2`$ or $`\tau = 0.1`$ it certifies at most 1% for the SAR-only generators. The cap does not help here, because the range of $`L - \epsilon K`$ is two-sided.

---

## 12. Computing everything from histograms

The GPU pass stores only sufficient statistics. For every patch and each of 200 score bins $`j`$ it keeps the number of pixels $`C_j`$ and of wrong pixels $`B_j`$, as cumulative counts. Then, for any threshold bin $`j`$:

```math
K(j) = \frac{C_j}{P}, \qquad \ell(j) = \frac{B_j}{P}.
```

**Capped loss.** For a cap $`q`$, let $`v`$ be the cumulative count of wrong pixels at the point where $`qP`$ pixels are trusted, linearly interpolated inside the bin that contains that point. Then

```math
\ell_q(j) = \frac{\min\{ B_j, v \}}{P}, \qquad K_q(j) = \frac{\min\{ C_j, qP \}}{P}.
```

This equals the expected loss under random tie-breaking inside that bin, so calibrating on these values certifies the randomized mask. All certificate experiments therefore run on a CPU in under a minute.

---

## 13. Experimental design

**Data.**

| Dataset | Use | Unit |
|---|---|---|
| SEN12MS-CR [5] (Sentinel-1 VV/VH, Sentinel-2) | generators, scores, calibration, winter hold-out | scene = one region in one season |
| SAR2Opt [6] (TerraSAR-X, 1 m) | second sensor | city (10 cities) |
| DFC2020 | external shift and land-cover check (supplement only) | none (no scene ids) |

**Preprocessing.** VV is clipped to $`[-25, 0]`$ dB and VH to $`[-32.5, -5]`$ dB. The B4, B3, B2 digital numbers are scaled by $`1/3000`$ into $`[0, 1]`$, so $`\tau = 0.1`$ equals 0.03 in reflectance.

**Scene-disjoint split.** 62 scenes (7,766 patches) train the generators; 22 train the trust scores; 64 scenes (64 patches each) calibrate and test; all 27 winter scenes are held out.

**Generators (all frozen before certification).** U-Net with an $`L_1`$ loss (54.4 M parameters); Pix2Pix, also with Monte Carlo dropout and with SAR-TopoGuard post-processing [7]; a rectified flow [8] (8.2 M parameters, 20 steps), as one sample and as the mean of five; the $`\tau`$-mode of the mixture density network; and a cloud-removal U-Net that also receives the cloudy optical image.

**Trust scores.** A learned verifier (13.6 M parameters) that predicts $`\mathbb{P}(\text{wrong})`$ from the SAR patch and the output; the floor; the spread of samples; and an edge-support score.

**Operating points.** Main: $`\alpha = 0.05`$, $`\tau = 0.1`$, $`\delta = 0.1`$. Secondary: $`\alpha = 0.1`$, $`\tau = 0.2`$.

**Two evaluations.**
- *Held-out splits:* 300 random splits of the 64 scenes into 32 for calibration and 32 for testing.
- *Scene superpopulation:* the 64 real scenes are treated as the population; $`m`$ scenes are drawn with replacement (500 draws), so the true risk $`R(\lambda)`$ of every threshold is known exactly.

**Pre-registration.** Hypotheses and pass criteria were written before each run. The notebook prints every verdict: 21 hypotheses in total, with seven failures reported honestly.

---

## 14. Results

**Table 1. SEN12MS-CR at $`\alpha = 0.05`$, $`\tau = 0.1`$** (32 calibration and 32 test scenes, 300 splits). "Over" is the share of splits whose test risk exceeds $`\alpha`$.

| Generator | PSNR (dB) | Patch CRC over | Patch CRC risk | Scene CRC over | Scene CRC area | Exact bound, uncapped, area | Capped certificate over | Capped certificate area | Error among trusted |
|---|---|---|---|---|---|---|---|---|---|
| $`\tau`$-mode (density net) | 19.98 | 0.51 | 0.055 | 0.20 | 0.097 | 0.000 | 0.04 | 0.086 | 0.29 |
| U-Net ($`L_1`$) | 19.78 | 0.45 | 0.050 | 0.02 | 0.119 | 0.000 | 0.03 | 0.098 | 0.31 |
| Pix2Pix | 19.67 | 0.49 | 0.053 | 0.07 | 0.127 | 0.000 | 0.05 | 0.102 | 0.27 |
| Pix2Pix, MC dropout | 19.67 | 0.48 | 0.053 | 0.06 | 0.128 | 0.000 | 0.05 | 0.102 | 0.27 |
| Pix2Pix + SAR-TopoGuard | 19.72 | 0.48 | 0.053 | 0.06 | 0.129 | 0.000 | 0.05 | 0.103 | 0.27 |
| Rectified flow, sample | 18.31 | 0.41 | 0.050 | 0.01 | 0.066 | 0.000 | 0.02 | 0.085 | 0.34 |
| Rectified flow, mean | 19.10 | 0.46 | 0.049 | 0.00 | 0.057 | 0.000 | 0.03 | 0.080 | 0.36 |
| Cloud removal U-Net (also sees the cloudy image) | 22.92 | 0.49 | 0.052 | 0.07 | 0.171 | 0.000 | 0.05 | 0.141 | 0.20 |

**Capped certificates.** Over 96 superpopulation configurations (eight generators, three operating points, $`m \in \{8, 16, 32, 64\}`$), the worst share of users over budget is 2.6% (target 10%), and the bound of Proposition 3 is never violated. For Pix2Pix, the capped certificate certifies 7.4% of pixels from 8 scenes and 8.9% from 16; the uncapped exact bound certifies none below 64 scenes.

**Confirmation on held-out winter scenes** (same rule and constants, no tuning). In 72 configurations the worst share over budget is 0.8%. The capped certificate always certifies more than the trivial area $`\alpha`$ (6.7%–34.8%); the uncapped bound never does.

**Table 2. Few labelled scenes.** Pix2Pix at $`\alpha = 0.05`$, $`\tau = 0.1`$ on winter ($`k`$ labelled scenes, 300 splits), and SAR2Opt at $`\alpha = 0.1`$ (five calibration and five test cities, 252 splits).

| Setting | Uncapped CRC area | Capped CRC risk | Capped CRC area | Capped certificate over | Capped certificate area |
|---|---|---|---|---|---|
| Winter, $`k = 4`$ | 0.000 | 0.038 | 0.083 | 0.00 | 0.058 |
| Winter, $`k = 8`$ | 0.000 | 0.037 | 0.094 | 0.01 | 0.066 |
| Winter, $`k = 12`$ | 0.000 | 0.038 | 0.104 | 0.03 | 0.073 |
| Winter, $`k = 16`$ | 0.000 | 0.039 | 0.112 | 0.03 | 0.078 |
| SAR2Opt, edge score | 0.000 | 0.074 | 0.110 | 0.00 | 0.088 |
| SAR2Opt, dropout spread | 0.000 | 0.049 | 0.071 | 0.00 | 0.060 |

**Ideas that did not work** (all re-run in Part IV of the notebook):
- Using the verifier's own predicted loss on unlabelled scenes as a control variate: it correlates only 0.29–0.67 with the true scene loss.
- Trusting the best share $`q`$ of *every* patch: the best 1% of pixels in a patch are wrong 41% of the time, against 50% for all pixels. The scores rank patches far better than they rank pixels inside a patch.
- A skew-corrected Student-$`t`$ bound: 8%–28% of users over budget at $`\delta = 0.1`$.
- Choosing the cap from a grid with a Bonferroni correction: valid, but it certified less than the fixed rule in 25 of 27 configurations.
- Bets $`c = 0.5`$ and $`c = 1`$: $`c = 0.5`$ certifies up to 0.07 less; $`c = 1`$ is within 0.02 of $`c = 0.75`$ but over budget up to 3% of the time.

---

## 15. What is certified, and what is not

- **The guarantee is on false trust as a share of *all* pixels.** Among trusted pixels, 20%–36% are still wrong at the strict tolerance $`\tau = 0.1`$ (Table 1, last column). Trusting any $`\alpha`$-share of pixels is valid without labels, so every table reports the certified area next to $`\alpha`$.
- **Scenes must be exchangeable.** A new region or season that differs from the calibration scenes needs new labelled scenes; the winter results show that 4–8 of them suffice with a cap.
- **Concept shift cannot be certified without labels.** If the optical appearance changes while the radar does not, no label-free method can tell.
- **The cap constants are exploratory.** They were chosen on the calibration scenes; only the winter results are confirmatory.

---

## 16. Running the notebook

The notebook [`CERTEO.ipynb`](CERTEO.ipynb) runs the whole study on real data only. It is written for [molab](https://molab.marimo.io) (marimo's hosted notebooks) with one GPU, and runs in about 2.5 hours on one RTX PRO 6000. It also runs in Jupyter.

**Parts.**

| Part | Cells | What it does |
|---|---|---|
| I | 1–14 | data access, models, training, evaluation, pre-registered hypotheses H1–H7, smoke and full runs |
| II | 15–24 | theory checks (H8–H13) and a dump of every result |
| III | 25–37 | certification tax, capped masks, product certificate, winter and SAR2Opt confirmation (H14–H20) |
| IV | 38–41 | ideas that did not work, qualitative example, both paper figures |

**Storage.** By default every result is written to `./certeo_work`. To keep results across sessions in your own Hugging Face bucket, set these environment variables or notebook secrets (no credential is stored in the code):

| Variable | Meaning |
|---|---|
| `HF_TOKEN` | Hugging Face token; only needed for a private results bucket |
| `CERTEO_BUCKET` | `<owner>/<bucket>` of your own Hugging Face bucket (optional) |
| `CERTEO_RUN_ROOT` | folder inside the bucket (default `certeo_v1`) |
| `CERTEO_OFFLINE=1` | force local-only storage |

**Resume.** Every stage writes atomic checkpoints. If a session ends, run all cells again: finished work is found and skipped.

**Resources.** Peak memory is about 15.4 GB of RAM and 8.0 GB of GPU memory. Datasets are streamed from their public Hugging Face mirrors and are never loaded whole.

**Outputs.** `HYPOTHESIS_REPORT.md`, `PART3_REPORT.md`, `PART4_REPORT.md`, JSON files for every table, and the figures `paper/fig_curves.pdf` and `paper/fig_qual.pdf`.

---

## 17. Repository layout

```
CERTEO.ipynb            the complete study (43 cells, Parts I–IV)
docs/certeo_method.tex  the mathematics of this README as a LaTeX document
docs/certeo_method.pdf  compiled version
README.md               this file
.gitignore              keeps local results and LaTeX build files out of git
LICENSE                 MIT
```

---

## 18. References

1. A. N. Angelopoulos, S. Bates, A. Fisch, L. Lei, and T. Schuster, "Conformal risk control," *ICLR*, 2024.
2. P. Hall, *The Bootstrap and Edgeworth Expansion*. Springer, 1992.
3. I. Waudby-Smith and A. Ramdas, "Estimating means of bounded random variables by betting," *J. R. Stat. Soc. B*, 2024.
4. C. K. Chow, "On optimum recognition error and reject tradeoff," *IEEE Trans. Inf. Theory*, 1970.
5. P. Ebel, A. Meraner, M. Schmitt, and X. X. Zhu, "Multisensor data fusion for cloud removal in global and all-season Sentinel-2 imagery," *IEEE Trans. Geosci. Remote Sens.*, 2021.
6. Y. Zhao, T. Celik, N. Liu, and H.-C. Li, "A comparative analysis of GAN-based methods for SAR-to-optical image translation," *IEEE Geosci. Remote Sens. Lett.*, 2022.
7. SAR-TopoGuard: training-free, speckle-aware topological verification for structure-preserving SAR-to-optical translation. *IEEE Geosci. Remote Sens. Lett.*, to appear.
8. X. Liu, C. Gong, and Q. Liu, "Flow straight and fast: learning to generate and transfer data with rectified flow," *ICLR*, 2023.
9. R. R. Bahadur and L. J. Savage, "The nonexistence of certain statistical procedures in nonparametric problems," *Ann. Math. Statist.*, 1956.
10. R. Dunn, L. Wasserman, and A. Ramdas, "Distribution-free prediction sets for two-layer hierarchical models," *J. Amer. Statist. Assoc.*, 2023.

## License

MIT; see [`LICENSE`](LICENSE). The datasets keep their own licences (SEN12MS-CR: CC BY 4.0; SAR2Opt: MIT).
