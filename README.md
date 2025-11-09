<a id="top"></a>

#  Contents

- [1. Choosing the Time Scale and Handling Delayed Entry in Cox Models](#sec-time-scale-delayed-entry)
- [2. Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model](#sec-adjusted-survival)
- [3. Survival Prediction Using Machine Learning](#survival-prediction-using-machine-learning)

---

<a id="sec-time-scale-delayed-entry"></a>
# 1. Choosing the Time Scale and Handling Delayed Entry in Cox Models


## Data Structure


### 📆 Calendar Timeline

- **Registry start:** **2000**  
  Cancer incidence recording begins.

- **Cohort entry window:** **2010–2030**  
  The analytic entry window (2010–2030) defines when individuals become eligible for inclusion—  
  that is, when covariates are collected, or when a person first becomes observable for follow-up.

- **Follow-up period:** **Entry date → 2050**  
  Each participant is followed from their entry date (which may occur years after diagnosis) until death, loss to follow-up, or study end.

- **Study end:** **2050**  
  Individuals alive without the event by this date are censored at their 2050 status.

---

### 👥 Inclusion Criteria

- Individuals **alive and observable at any point between 2010 and 2030** are included,  
  regardless of diagnosis year, as long as they survived long enough to enter the analytic window.

> **Note:** This cohort design uses **delayed entry (left truncation)** — follow-up begins only once a participant becomes observable (their entry date), rather than at the registry start in 2000.  
> This approach ensures that survival time is counted only when participants are under observation, avoiding **immortal time bias** — an artificial inflation of survival that occurs when unobserved years before entry are mistakenly included.
>
> **Example — Immortal Time Bias**
>
> A person is **diagnosed in 2008**, but your study begins observation in **2010**.  
> If survival is measured from 2008, that person automatically appears to have survived at least **2 years**, because anyone who died before 2010 would never be observed.  
> The period **2008–2010** becomes “**immortal time**” — time during which the person was guaranteed to survive simply to be included in the study.


### Outcomes and Censoring
- **Event of interest:** death from any cause (overall survival).  
- **Censoring:** alive at 2050 or lost to follow-up.  
- **Status variable:** `status = 1` for death from any cause, `status = 0` otherwise.

## ⏳ Choosing the Time Scale

In survival analysis, the time scale defines how follow-up is measured.
Two standard choices are time since diagnosis and attained age.

## Option 1. Time Since Diagnosis (Time-on-Study)

### Definition:
Follow-up starts at the diagnosis date (t = 0), and time measures years survived after diagnosis.

### Use when:

The focus is on survival after diagnosis (e.g., 5-year survival).

You want time measured relative to diagnosis, not absolute age.

### Setup:
```{r}
start_time = max(0, years_between(entry_date, diag_date))
stop_time  = years_between(event_date, diag_date)
```

### Interpretation:
S(2) = probability of surviving beyond 2 years after diagnosis.
age_at_dx (age at diagnosis) can be included as a baseline covariate.

### Example model:

```{r}
fit <- coxph(Surv(t_start, t_stop, status) ~ age_at_dx + gene + other_covs, data = d)
survfit(fit)
```

## Option 2. Attained Age (Age as Time Scale)

### Definition:
The analysis clock measures chronological age,
where attained_age = age_at_diagnosis + time_since_diagnosis.

### Use when:

The interest is in how risk changes with age itself.

You want survival expressed directly by age (e.g., “probability of surviving beyond age 70”).

### Setup:
```{r}
start_age = max(age_at_entry, age_at_diagnosis)
stop_age  = age_at_event
```

### Interpretation:
The survival function S(age) gives the probability of remaining event-free up to and beyond a given age,
conditional on being event-free at entry.

### Example:

S(40) = 0.85 → among participants who were event-free at entry and share the same covariates,
85% are expected to survive beyond age 40.

---

## References

* Cox DR. 1972. Regression Models and Life-Tables. JRSS B 34:187–220.
* Therneau TM, Grambsch PM. 2000. Modeling Survival Data: Extending the Cox Model. Springer.
* Klein JP, Moeschberger ML. 2003. Survival Analysis: Techniques for Censored and Truncated Data. Springer.


---

<a id="sec-adjusted-survival"></a>
# 2. Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model


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
# 3. Survival Prediction Using Machine Learning

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

