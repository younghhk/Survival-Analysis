# Adjusted Survival Curve Estimation from the Cox Proportional Hazards Model

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


