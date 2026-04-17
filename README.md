# Biological_time_series_analysis

Introduction to Time-Series Analysis of Birth Dynamics
This report presents a comprehensive time-series analysis of birth counts—indicative of infestation levels—across two distinct macro-areas: Susegana and Breganze.

Key Analytical Insights:Statistical Foundation: Given the strict nature of count data within a fixed time interval, the target variable was formally modeled utilizing a Poisson distribution.

Temporal Correlations: The analysis reveals a strong autoregressive correlation in both locations, indicating that the state of infestation at a given time is highly predictive of the subsequent state. Susegana exhibits a uniform, easily identifiable epidemic peak. Conversely, Breganze displays asynchronous dynamics and a distinct rebound effect around day 60.

Methodological Advancements: Classical polynomial models consistently underestimate the maximum epidemiological peaks, which are heavily concentrated between days 15 and 25. To resolve this, a Generalized Linear Mixed Model (GLMM) was employed.

Primary Findings: The advanced GLMM successfully incorporates an AR1 correlation structure to capture biological memory. 
Furthermore, it demonstrates that the birth rate in Breganze is significantly lower, representing approximately 8.6% of Susegana's rate. Ultimately, this framework provides superior management of overdispersion and temporal dependencies compared to standard parametric models.
