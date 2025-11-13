<a id="top"></a>



[![Cancer Research Software Hub](https://img.shields.io/badge/Back_to-Hub-blue)](https://github.com/younghhk/NCI)



# 📘 Contents


- [Handling Delayed Entry in Cox Models](#sec-time-scale-delayed-entry)
- [Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model](#sec-adjusted-survival)
- [Survival Prediction Using Machine Learning](#survival-prediction-using-machine-learning)

---

<a id="sec-time-scale-delayed-entry"></a>
#  Handling **Delayed Entry** (Left Truncation) in Cox Models

---

## 1) Why Left Truncation (Delayed Entry) Must Be Handled Explicitly

In survival analysis, each subject *i* contributes information to the likelihood **only** for the period during which they are **actually under observation**.  
If some participants enter the study late (e.g., due to recruitment between 2010–2030), their earlier, unobserved time cannot contribute valid risk time — otherwise, survival will be overestimated (**immortal time bias**).

### Conceptual Framework

Let:
- **Tᵢ** — true event time for subject *i*  
- **Lᵢ** — entry time (left truncation point)  
- **Cᵢ** — censoring time  
- Subject *i* is **observable only if** Tᵢ > Lᵢ  
- Observed data: (Lᵢ, Tᵢ, δᵢ, Xᵢ), where   δᵢ = 1 if the event is observed, and 0 otherwise.
  



In the Cox model, each time an event occurs for subject *i*, we compare that subject’s covariates to all others who were **at risk** at that same moment.

The **risk set** at the event time of subject *i* is:

> **R(Tᵢ) = { j : Lⱼ < Tᵢ ≤ Tⱼ }**

where:
- **i** indexes the subject who experienced the event at time **Tᵢ**,  
- **j** indexes *all* subjects who were still “under observation” at that moment —  
  meaning they had already entered the study (**Lⱼ < Tᵢ**) and had not yet experienced the event or censoring (**Tⱼ ≥ Tᵢ**).  
- **δᵢ = 1** if subject *i* experienced the event, and **0** if they were censored.



The **Cox proportional hazards model** assumes:

> λ(t | Xᵢ) = λ₀(t) × exp(Xᵢᵀβ)

The **partial likelihood** is then:

> L(β) = ∏ᵢ [ exp(Xᵢᵀβ) / Σⱼ∈R(Tᵢ) exp(Xⱼᵀβ) ]^δᵢ

Each event time contributes one term:  
- The **numerator** corresponds to the subject who actually experienced the event (*i*).  
- The **denominator** sums the risk contributions of all subjects *j* who were **at risk** at that moment.


If delayed entry is ignored, the risk set is incorrectly defined as:

> **R(Tᵢ) = { j : Tⱼ ≥ Tᵢ }**

This version includes people who **had not yet entered** the study by time **Tᵢ** (i.e., Lⱼ ≥ Tᵢ),  
artificially inflating the denominator and biasing hazard estimates.

By correctly using **R(Tᵢ) = { j : Lⱼ < Tᵢ ≤ Tⱼ }**,  
we ensure that subjects contribute **only during the period they are actually observable**.  



###  Model Implementation in R

In `R`, left truncation is implemented using the `(start, stop]` form of `Surv()`:

```r
coxph(Surv(start_time, stop_time, status) ~ covariates, data = d)
```

---
## 2) Study Frame & Data Structure (Hypothetical)

### 📆 Calendar Timeline
- **Registry start:** **2000** – cancer cases begin being recorded.
- **Cohort entry window:** **2010–2030** – people become **observable/eligible** (e.g., recruited, baseline covariates recorded, or first linked).
- **Follow-up:** **Entry date → 2050** – death or censoring observed.


### 👥 Inclusion (Who gets in?)
- Individuals **alive and observable at any time between 2010 and 2030**, regardless of diagnosis year, **provided they survived to be observed**.
  
> **Why delayed entry matters (a.k.a. left truncation):**  
> If someone was diagnosed in **2008** but you only begin observing them in **2012**, they must have survived until 2012 to appear in your data.  
> Counting **2008–2012** as time-at-risk would add **“immortal time”** (time they had to survive to be included), biasing survival **upward**.  
> **Solution:** Start their risk time at the **entry date** (2012), not the registry start.

### Outcomes & Censoring (for this demo)
- **Event of interest:** **death from any cause** (overall survival).
- **Censoring:** alive at 2050 or lost to follow-up.
- **Status variable:** `status = 1` (death), `status = 0` (censored).

---

## 3) Two Valid Time Scales

### Option A. **Time Since Diagnosis** (time-on-study)
- **Definition:** Clock starts at **diagnosis** (`t = 0`), time measures **years after diagnosis**.
- **Use when:** You want “**5-year survival** after diagnosis” or durations relative to diagnosis.
- **Left truncation:** risk starts at the **later** of (0, time from diagnosis to entry).

### Option B. **Attained Age** (age as time scale)
- **Definition:** Clock is **chronological age**; `attained_age = age_at_diagnosis + time_since_diagnosis`.
- **Use when:** Risk is primarily **age-driven** and you want survival **by age** (e.g., `S(70)`).
- **Left truncation:** risk starts at the **later** of (age at entry, age at diagnosis).

---

## 4) Construct (start, stop] times with delayed entry

- Time-since-diagnosis: t_start = max(0, entry - diagnosis), t_stop = event - diagnosis.

- Attained age: age_start = max(age_at_entry, age_at_dx), age_stop = age_at_event.

## 5) Fit Cox Models With Left Truncation

### A) Time Since Diagnosis (time-on-study)


```{r}
fit_dx <- coxph(Surv(t_start, t_stop, status) ~ age_at_dx + X, data = d)
summary(fit_dx)
```

**Interpretation**

- `S(2)` from this model = probability of surviving >2 years after diagnosis for the specified covariate profile.
- `age_at_dx` is a baseline covariate here (okay to include).

### B) Attained Age (age as time scale)

 When age is the time scale, the `Surv()` uses attained age at entry/exit
 
```{r}
fit_age <- coxph(Surv(age_start, age_stop, status) ~ X, data = d)
summary(fit_age)
```
**Interpretation**

- `S(70)` = probability of being event-free at age 70, conditional on being event-free at one’s entry age.

- Do not add a generic “age” covariate when age is the time axis.



---

<a id="sec-adjusted-survival"></a>
# Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model


This repository provides a reproducible example of how to generate **adjusted survival curves** from a **Cox proportional hazards model** in R. These curves resemble Kaplan–Meier plots but are model-based estimates that account for covariates. They are useful for visualizing survival probabilities while controlling for confounding factors.



###  Setup

Install R (≥ 4.0.0) and the following packages:

```r
install.packages(c("survival", "survminer"))
````



### Data Example

Your dataset should contain at least:

* `time`: follow-up or survival time
* `status`: event indicator (1 = event, 0 = censored)
* Covariates: explanatory variables (e.g., treatment group, age)

Example dataset structure:

```r
head(mydata)
#   time status treatment age
# 1   10      1        A  65
# 2   15      0        B  72
# 3   20      1        A  54
```



###  Code Example

#### 1. Fit the Cox Model

```r
library(survival)
library(survminer)

# Fit Cox proportional hazards model
cox_fit <- coxph(Surv(time, status) ~ treatment + age, data = mydata)
summary(cox_fit)
```


#### 2. Generate Adjusted Survival Estimates

```r
# Survival curve for treatment A, age 65
cox_surv <- survfit(cox_fit, newdata = data.frame(treatment = "A", age = 65))

# Plot adjusted survival curve
ggsurvplot(cox_surv, 
           conf.int = TRUE, 
           ggtheme = theme_minimal(),
           title = "Adjusted Survival Curve (Cox Model)")
```



#### 3. Compare Two Groups

```r
# Survival estimates for treatment A vs B, both at age 65
newdata_list <- list(
  data.frame(treatment = "A", age = 65),
  data.frame(treatment = "B", age = 65)
)

cox_surv_groups <- survfit(cox_fit, newdata = newdata_list)

ggsurvplot(cox_surv_groups, 
           conf.int = TRUE, 
           ggtheme = theme_minimal(),
           legend.title = "Treatment Group",
           title = "Cox-Adjusted Survival Curves")
```



###  R Example



```r
# Load dataset
data(lung)

# Fit Cox model with sex + age
cox_fit <- coxph(Surv(time, status) ~ sex + age, data = lung)

# Adjusted survival curve for a 60-year-old male
ggsurvplot(survfit(cox_fit, newdata = data.frame(sex = 1, age = 60)),
           conf.int = TRUE,
           ggtheme = theme_minimal(),
           title = "Adjusted Survival Curve (Cox Model, lung dataset)")
```



###  Interpretation

* **Curve meaning**: Probability of surviving past a given time, adjusted for covariates in the model.
* **Confidence intervals**: Shaded bands represent estimation uncertainty.
* **Group comparison**: Separation between curves shows survival differences after covariate adjustment.
* **Hazard ratios**: From `summary(cox_fit)`. HR > 1 = higher risk, HR < 1 = lower risk compared to reference.

⚠️ These curves are **model-based estimates**, not raw Kaplan–Meier curves. They provide a clearer picture of covariate-adjusted survival differences.



###  References

* Cox DR (1972). *Regression Models and Life-Tables*. Journal of the Royal Statistical Society, Series B, 34(2):187–220.
* Therneau TM (2023). **A Package for Survival Analysis in R**. [CRAN survival package](https://cran.r-project.org/package=survival)
* Kassambara A, Kosinski M, Biecek P (2023). **survminer: Drawing Survival Curves using 'ggplot2'**. [CRAN survminer package](https://cran.r-project.org/package=survminer)

---
<a id="survival-prediction-using-machine-learning"></a>
# Survival Prediction Using Machine Learning

###  Examples with allogeneic hematopoietic cell transplantation (HCT):

- **Model Performance (AUC Comparison):** [AUC.pdf](model_auc_comparison.pdf)  
- **Survival Curves (Kaplan–Meier):** [KM.pdf](km_plots.pdf)  
- **Risk Stratification (Random Survival Forests):** [Risk Stratification.pdf](rsf_risk_groups.pdf)  

###  Resources
- [Random Survival Forests Overview](https://www.randomforestsrc.org/articles/survival.html)  

### 🔒 Repository Access

The **Survival Prediction Using Machine Learning** materials use examples derived from a restricted dataset and are available to **approved collaborators only**.

Access requests may be directed to:

**Grace Hong**  
📧 [grace.hong@nih.gov](mailto:grace.hong@nih.gov)

*(Please include your GitHub username and affiliation in your request.)*

