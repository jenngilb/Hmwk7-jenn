import SurvivalModelClasses as Cls
import scr.SamplePathClasses as SamplePathSupport
import scr.FigureSupport as Fig
import scipy.stats as stat

MORTALITY_PROB = 0.1    # annual probability of mortality
TIME_STEPS = 100        # simulation length
SIM_POP_SIZE = 1000     # population size of the simulated cohort
ALPHA = 0.05            # significance level

# create a cohort of patients
myCohort = Cls.Cohort(id=1, pop_size=SIM_POP_SIZE, mortality_prob=MORTALITY_PROB)

# simulate the cohort
cohortOutcome = myCohort.simulate(TIME_STEPS)

# plot the sample path
SamplePathSupport.graph_sample_path(
    sample_path=cohortOutcome.get_survival_curve(),
    title='Survival Curve',
    x_label='Time-Step (Year)',
    y_label='Number Survived')

# plot the histogram
Fig.graph_histogram(
    data=myCohort.get_survival_times(),
    title='Histogram of Patient Survival Time',
    x_label='Survival Time (Year)',
    y_label='Count')

#Hmwk Problem 1
print("Hmwk Q1: The cohort survival probability past 5 years is", myCohort.cohortsurvivalresults())

#Hmwk Problem 2
print("Hmwk Q2: The type of distribution is binomial, and the only parameter is the patient mortality.")

#Hmwk Problem 3
print("Hmwk Q3:The likelihood ratio is:", stat.binom.pmf(k=400,n=573,p=0.5))

#Hmwk Problem 4

MORTALITY_PROB = 573/1000    # annual probability of mortality
TIME_STEPS = 5        # simulation length
REAL_POP_SIZE = 573     # size of the real cohort to make the projections for
NUM_SIM_COHORTS = 1000   # number of simulated cohorts used for making projections
ALPHA = 0.05            # significance level

# calculating prediction interval for mean survival time
# create multiple cohorts
multiCohort = Cls.MultiCohort(
    ids=range(NUM_SIM_COHORTS),   # [0, 1, 2 ..., NUM_SIM_COHORTS-1]
    pop_sizes=[REAL_POP_SIZE] * NUM_SIM_COHORTS,  # [REAL_POP_SIZE, REAL_POP_SIZE, ..., REAL_POP_SIZE]
    mortality_probs=[MORTALITY_PROB]*NUM_SIM_COHORTS  # [p, p, ....]
)
# simulate all cohorts
multiCohort.simulate(TIME_STEPS)

# print projected mean survival time (years)
print('Hmwk Q4: Projected mean survival time (years) is',
      multiCohort.get_overall_mean_survival())
# print projection interval
print('95% projection interval of average survival time (years)',
      multiCohort.get_PI_mean_survival(ALPHA))

# print the patient survival time
print('Average survival time (years):', cohortOutcome.get_ave_survival_time())
print('95% CI of average survival time (years)', cohortOutcome.get_CI_survival_time(ALPHA))
# report mean and projection interval

#Hmwk Problem 5
print('Hmwk Q5: Mean survival time and {:.{prec}%} projection interval:'.format(1 - ALPHA, prec=0),
      Cls.calibrated_model.get_mean_survival_time_proj_interval(ALPHA, deci=4))

#Hmwk Problem 6

MORTALITY_PROB = 800/1146    # annual probability of mortality
TIME_STEPS = 5        # simulation length
REAL_POP_SIZE = 1146     # size of the real cohort to make the projections for
NUM_SIM_COHORTS = 1000   # number of simulated cohorts used for making projections
ALPHA = 0.05            # significance level

# calculating prediction interval for mean survival time
# create multiple cohorts
multiCohort = Cls.MultiCohort(
    ids=range(NUM_SIM_COHORTS),   # [0, 1, 2 ..., NUM_SIM_COHORTS-1]
    pop_sizes=[REAL_POP_SIZE] * NUM_SIM_COHORTS,  # [REAL_POP_SIZE, REAL_POP_SIZE, ..., REAL_POP_SIZE]
    mortality_probs=[MORTALITY_PROB]*NUM_SIM_COHORTS  # [p, p, ....]
)
# simulate all cohorts
multiCohort.simulate(TIME_STEPS)

# print projected mean survival time (years)
print('Hmwk Q6: Projected mean survival time (years) is',
      multiCohort.get_overall_mean_survival())
# print projection interval
print('95% projection interval of average survival time (years)',
      multiCohort.get_PI_mean_survival(ALPHA))
# Estimate of mortality probability and the posterior interval
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1-ALPHA, prec=0),
      Cls.get_mortality_estimate_credible_interval(ALPHA, 4))

print("My code did not work, but if it had I suspect the intervals for the transient state simulation would both have a smaller range because there would be a higher number of observations, as 1000>573. This is related to the Law of Large Numbers --> more observations can provide greater certainty than could otherwise be ascertained.")
print("Steady state models are used when the Law of Large numbers can be applied already, so my guess is that the projection interval WOULD have a more narrow range of values.")



