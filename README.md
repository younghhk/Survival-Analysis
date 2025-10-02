<a id="top"></a>

# 📚 Contents

- [1. Choosing the Time Scale and Handling Delayed Entry in Cox Models](#sec-time-scale-delayed-entry)
- [2. Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model](#sec-adjusted-survival)

---

<a id="sec-time-scale-delayed-entry"></a>
# 1. Choosing the Time Scale and Handling Delayed Entry in Cox Models


## Study Setup and Data Structure

### Calendar timeline

- **Registry start:** 2000  
  Incident cases begin to be recorded.

- **Cohort entry window:** **2010–2030**  
  Individuals are eligible to **enter** the analytic cohort on their first qualifying date within this window (e.g., recruitment/sampling date when covariates are recorded).

- **Follow-up period:** entry date (between 2010–2030) through **2050**  
  Events and censoring are recorded through 2050.

- **Study end:** 2050  
  Anyone without the event by this date is censored at their status in 2050.

### Who is included

- People who are alive and eligible **at their entry date between 2010 and 2030**, plus individuals diagnosed after 2010 who enter when they are first observed within the entry window.  
- People who were diagnosed before 2010 but **died before their potential entry date** (i.e., before 2010) are not in the analytic cohort.  
- Individuals diagnosed **after 2030** (and thus not observed within the 2010–2030 entry window) are **not** part of the cohort.

This is handled with **delayed entry (left truncation)** so that risk time begins only when a person is actually observable (their entry date) and, for post-diagnosis outcomes, after diagnosis.

### Outcomes and censoring

- **Event of interest:** disease-specific death (replace with your exact endpoint if different).  
- **Censoring:** alive at 2050, competing-cause death, or lost to follow-up.  
- **Status variable:** `status = 1` for the event of interest, `status = 0` otherwise.

### Two valid analysis choices for the time axis

1) **Time since diagnosis**  
   - The analysis time measures years since each person’s diagnosis.  
   - Useful for questions like “survival *t* years after diagnosis.”  
   - To report survival by *attained age*, you must also specify age at diagnosis and convert: `age = age_at_dx + t`.

2) **Attained age**  
   - The analysis time is chronological age.  
   - Useful for questions like “survival at a given age (20, 30, 40…).”  
   - Values like `S(30)` can be read directly at age 30.

Both choices require **delayed entry between 2010 and 2030** so each subject contributes person-time only after they actually enter observation, and (for post-diagnosis outcomes) only after diagnosis has occurred.



### How delayed entry works in this setup

A person contributes risk only after both conditions hold:
1) the person has been **diagnosed**, and  
2) the person is **observable** within the analytic **entry window (2010–2030)**.

- On the **time-since-diagnosis** scale:  
  `start_time = max(0, years_between(entry_date, diag_date))`  
  `stop_time  = years_between(event_date, diag_date)`

- On the **attained-age** scale:  
  `start_age = max(age_at_entry, age_at_dx)`  
  `stop_age  = age_at_event`

### Examples

- **Diagnosed 2008 at age 55; enters 2012 (age 59); dies 2016 (age 63)**  
  - Time since diagnosis: `start = 4` (2012−2008), `stop = 8` (2016−2008), `status = 1`  
  - Attained age: `start = 59`, `stop = 63`, `status = 1`

- **Diagnosed 2029 at age 62; enters 2029 (age 62); alive in 2050 (age 83)**  
  - Time since diagnosis: `start = 0`, `stop = 21`, `status = 0`  
  - Attained age: `start = 62`, `stop = 83`, `status = 0`

- **Diagnosed 1999 at age 45; died 2009 at age 55**  
  - Not in cohort: death occurred before any possible entry in 2010–2030.

- **Diagnosed 2032 at age 50; first observed 2032**  
  - Not in cohort: diagnosis/observation occurs **after** the entry window (post-2030).

### Why this setup matters

- Using delayed entry aligned to **each subject’s entry date (2010–2030)** matches what was truly observed and prevents immortal time bias.  
- Choosing the time axis that matches your estimand makes interpretation straightforward:
  - **Time since diagnosis** answers “how long after diagnosis.”  
  - **Attained age** answers “at what age.”





## Case A. Time since diagnosis as the time scale

**What the clock measures**  
- t = 0 at each subject’s diagnosis.  
- Event time is years since diagnosis until death or censoring.

**Can age at diagnosis be a covariate here?**  
- Yes. Age at diagnosis (a0) is a baseline characteristic, not the time axis.  
- If age at diagnosis affects hazard, include it in the model or use strata.

**How to read survival**  
- `S(2)` means the probability of surviving more than 2 years after diagnosis.  
- Survival can differ by age at diagnosis if you include it in the model.

**Delayed entry in 2010 on the diagnosis scale**  
- Entry time since diagnosis = (calendar 2010) − (calendar diagnosis date), truncated at 0.  
- Example: diagnosed in 2008 → entry at t = 2. Diagnosed in 2012 → negative value, truncate to 0 (they enter at diagnosis).

**R template**

```r
library(survival)

# Variables (example):
# diag_date: calendar date of diagnosis
# entry_date: calendar date of analytic entry (commonly in 2010)
# event_date: calendar date of event or censor
# status: 1 = target-cause death, 0 = censor (alive or other-cause death)
# gene: fixed covariate
# age_at_dx: age at diagnosis (baseline covariate)

to_years <- function(days) as.numeric(days) / 365.25

d$t_start <- pmax(0, to_years(difftime(d$entry_date, d$diag_date, units = "days")))
d$t_stop  <- to_years(difftime(d$event_date, d$diag_date, units = "days"))

# Keep only subjects observed after entry
d <- subset(d, t_stop > t_start)

fit_A <- coxph(Surv(t_start, t_stop, status) ~ age_at_dx + gene + other_covs, data = d)
summary(fit_A)

# Model-based survival curves t years after diagnosis
# For reporting at attained age a for someone with age_at_dx = a0:
# evaluate S(t) at t = a - a0
sf_A <- survfit(fit_A, newdata = data.frame(age_at_dx = 30, gene = 1, other_covs = ...))
````

**Reporting “by age 20, 30, 40” in Case A**

* You must fix age at diagnosis.
* Example: if a0 = 30, then age 40 corresponds to t = 10. Read S(10).
* If you want to compare a0 = 20 vs 30 vs 40, plot separate curves by a0 or include a0 in the model and evaluate at chosen values.

---

## Case B. Cox Models With Attained Age as the Time Scale


**When to use attained age**

Use attained age as the time scale when your scientific question is about the probability of being event-free **at a given age**. Examples:
- Survival at age 20, 30, 40.
- Comparing survival by genotype while controlling for current age through the time scale.




**Delayed entry (left truncation) on the age scale**

Your analytic follow-up starts in 2010. Subjects contribute risk only **after** both of the following are true:
1) They have been diagnosed.
2) They are under observation beginning in 2010.

Therefore, the entry age for analysis is:
```r

start_age = max(age_at_2010, age_at_dx)

```r
The exit age is:
```r

stop_age = age_at_event

````
A subject is at risk for ages in the interval `[start_age, stop_age)`.

**Examples**
- Diagnosed in 2008 at age 55, age in 2010 is 57 → start_age = max(57, 55) = 57
- Diagnosed in 2015 at age 62, age in 2010 is 57 → start_age = max(57, 62) = 62
- Diagnosed in 1995 at age 40, age in 2010 is 55 → start_age = max(55, 40) = 55



**Model specification in R**

```r
library(survival)

# d must contain: age_at_dx, age_at_2010, age_at_event, status, gene, other covariates

d$start_age <- pmax(d$age_at_2010, d$age_at_dx)
d$stop_age  <- d$age_at_event

# Keep subjects with observed follow-up after entry
d <- subset(d, stop_age > start_age)

# Cox PH with attained age as the time scale
fit_age <- coxph(Surv(start_age, stop_age, status) ~ gene + other_covs, data = d)
summary(fit_age)

# Model-based survival curve by attained age for a chosen covariate profile
new_prof <- data.frame(gene = 1, other_covs = ...)
sf_age <- survfit(fit_age, newdata = new_prof)

# You can read survival at age 20, 30, 40, etc. directly from sf_age
````

Notes:

* Since attained age is the time scale, do **not** include a generic “age” covariate. It would duplicate the role of the baseline hazard.
* You may include **age at diagnosis** as a covariate if the scientific question is about earlier versus later diagnosis conditional on current age. Use with care because it is related to entry.



**How to interpret**

* The survival function (` S(a) `) represents the probability of being event-free **at age ( a )** for the specified covariate profile, **conditional on being event-free at the entry age** imposed by delayed entry.
* Hazard ratios from `summary(fit_age)` compare instantaneous risk at the same attained age and covariates.

Example interpretation:

* If (` S(40) = 0.85` ) for a given profile, the model estimates that 85% are event-free at age 40 among those who were event-free when they entered observation and who have that covariate profile.

---




## References

* Cox DR. 1972. Regression Models and Life-Tables. JRSS B 34:187–220.
* Therneau TM, Grambsch PM. 2000. Modeling Survival Data: Extending the Cox Model. Springer.
* Klein JP, Moeschberger ML. 2003. Survival Analysis: Techniques for Censored and Truncated Data. Springer.


---

<a id="sec-adjusted-survival"></a>
# 2. Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model


This repository provides a reproducible example of how to generate **adjusted survival curves** from a **Cox proportional hazards model** in R. These curves resemble Kaplan–Meier plots but are model-based estimates that account for covariates. They are useful for visualizing survival probabilities while controlling for confounding factors.



##  Setup

Install R (≥ 4.0.0) and the following packages:

```r
install.packages(c("survival", "survminer"))
````



## Data Example

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



##  Code Example

### 1. Fit the Cox Model

```r
library(survival)
library(survminer)

# Fit Cox proportional hazards model
cox_fit <- coxph(Surv(time, status) ~ treatment + age, data = mydata)
summary(cox_fit)
```


### 2. Generate Adjusted Survival Estimates

```r
# Survival curve for treatment A, age 65
cox_surv <- survfit(cox_fit, newdata = data.frame(treatment = "A", age = 65))

# Plot adjusted survival curve
ggsurvplot(cox_surv, 
           conf.int = TRUE, 
           ggtheme = theme_minimal(),
           title = "Adjusted Survival Curve (Cox Model)")
```



### 3. Compare Two Groups

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



##  Quick Start with Built-in Data

You can test the workflow with R’s built-in `lung` dataset:

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



##  Interpretation

* **Curve meaning**: Probability of surviving past a given time, adjusted for covariates in the model.
* **Confidence intervals**: Shaded bands represent estimation uncertainty.
* **Group comparison**: Separation between curves shows survival differences after covariate adjustment.
* **Hazard ratios**: From `summary(cox_fit)`. HR > 1 = higher risk, HR < 1 = lower risk compared to reference.

⚠️ These curves are **model-based estimates**, not raw Kaplan–Meier curves. They provide a clearer picture of covariate-adjusted survival differences.



##  References

* Cox DR (1972). *Regression Models and Life-Tables*. Journal of the Royal Statistical Society, Series B, 34(2):187–220.
* Therneau TM (2023). **A Package for Survival Analysis in R**. [CRAN survival package](https://cran.r-project.org/package=survival)
* Kassambara A, Kosinski M, Biecek P (2023). **survminer: Drawing Survival Curves using 'ggplot2'**. [CRAN survminer package](https://cran.r-project.org/package=survminer)


