
## EXPLANATION OF DATA EXAMPLES
# 

### SIMPLE TEST CASE 
# A small case with 3 objects and 3 judges to test whether all functions are
# working as intended.
# 

### SIMULATED DATA 1
# Rankings: 3x (1,2,3,4), 2x (1,2,4,3), 2x (2,1,3,4), 2x (2,1,4,3), 1x (3,2,4,1)
# Ratings: each unique
# The rankings generally represent equal support split between the first two 
# objects, which are almost always ranked as the best two. The last two objects
# are similarly equally supported but both recognized as significantly worse
# that the first two. 
# The ratings generally represent a consensus that the first object is better,
# but judges who don't rate it the best have stronger opinions against it. 
# There is even better consensus that the third object is better than the 
# fourth. 
# 

### SIMULATED DATA 2
# Stronger disagreement between rankings and ratings
# In particular, rankings indicate that 2 > 1 and 4 > 3, whereas ratings 
# indicate 1 > 2 and 3 > 4.
# 

### SALMON DATA FROM SADER ET AL.
# Data taken from Monday salmon samples only. 
# Rankings: (1,2,3,4) = (A3, B4, C5, D6)
# Likert scale is adjusted to ratings
#

## TESTING

# Load libraries, resources, and functions from file
source("conor/thesis/semester2/estimation-resources.R")

### SIMPLE TEST CASE
# Create ranking and rating data
rankings <- matrix(c(c(1,2,3),c(1,2,3),c(1,3,2)),nrow=3,byrow=T)
ratings <- matrix(c(c(1,1,1),c(0,1,2),c(0,1,2)),nrow=3,byrow=T)

# Select M and alpha
M <- 5
test_alpha <- 0.15

# Visualize the ranking and rating data 
test_rankings_plot <- plot_rankings(rankings)
test_ratings_plot <- plot_ratings(ratings)

# Calculate original MB and new WMB maximum likelihood estimators
test_MLE_mb <- show_mb_estimate(rankings, ratings, M)
test_MLE_wmb <- show_wmb_estimate(rankings, ratings, M, alpha=test_alpha)

# Create plots showing how estimation results change with alpha
test_p_plot <- ASTAR_plot(rankings,ratings,M,n=300,parameter="p",return_plot=T)
test_theta_plot <- ASTAR_plot(rankings,ratings,M,n=300,parameter="theta",return_plot=T)

# Save all plots to files
print("Saving test plot to conor/thesis/semester2/outputs/test_plot.png")
png("conor/thesis/semester2/outputs/test_rankings_plot.png", width=6, height=5, units="in", res=1200)
test_rankings_plot
dev.off()
png("conor/thesis/semester2/outputs/test_ratings_plot.png", width=6, height=5, units="in", res=1200)
test_ratings_plot
dev.off()
png("conor/thesis/semester2/outputs/test_p_plot.png", width=6, height=5, units="in", res=1200)
test_p_plot
dev.off()
png("conor/thesis/semester2/outputs/test_theta_plot.png", width=6, height=5, units="in", res=1200)
test_theta_plot
dev.off()


### SIMULATED DATA 1
# Create ranking and rating data
rankings <- matrix(c(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(1,2,4,3),c(1,2,4,3),c(2,1,3,4),c(2,1,3,4),c(2,1,4,3),c(2,1,4,3),c(3,2,4,1)),nrow=10,byrow=T)
ratings <- matrix(c(c(0,0,2,4),c(0,0,3,4),c(0,1,2,2),c(0,1,3,3),c(0,1,2,3),c(0,1,2,4),c(0,2,4,3),c(1,1,2,3),c(2,0,3,4),c(2,1,4,3)),nrow=10,byrow=T)

# Select M and alpha
M <- 4
sim1_alpha <- 0.15

# Visualize the ranking and rating data 
sim1_rankings_plot <- plot_rankings(rankings)
sim1_ratings_plot <- plot_ratings(ratings)

# Calculate original MB and new WMB maximum likelihood estimators
sim1_MLE_mb <- show_mb_estimate(rankings, ratings, M)
sim1_MLE_wmb <- show_wmb_estimate(rankings, ratings, M, alpha=sim1_alpha)

# Create plots of confidence intervals on p
sim1_CI_wmb <- ci_wmb(rankings=rankings,ratings=ratings,M=M,alpha=sim1_alpha,
                      interval=0.95,nsamples=200)
sim1_ci_plot <- plot_p_wmb(sim1_MLE_wmb, sim1_CI_wmb, show_data=T)

# Create plots showing how estimation results change with alpha
sim1_p_plot <- ASTAR_plot(rankings,ratings,M,n=100,parameter="p",return_plot=T)
sim1_theta_plot <- ASTAR_plot(rankings,ratings,M,n=100,parameter="theta",return_plot=T)

# Save all plots to files
png("conor/thesis/semester2/outputs/sim1_rankings_plot.png", width=6, height=5, units="in", res=1200)
sim1_rankings_plot
dev.off()
png("conor/thesis/semester2/outputs/sim1_ratings_plot.png", width=6, height=5, units="in", res=1200)
sim1_ratings_plot
dev.off()
png("conor/thesis/semester2/outputs/sim1_p_plot.png", width=6, height=5, units="in", res=1200)
sim1_p_plot
dev.off()
png("conor/thesis/semester2/outputs/sim1_theta_plot.png", width=6, height=5, units="in", res=1200)
sim1_theta_plot
dev.off()
png("conor/thesis/semester2/outputs/sim1_ci_plot.png", width=6, height=5, units="in", res=1200)
sim1_ci_plot
dev.off()


### SIMULATED DATA 2
rankings <- matrix(c(c(1,2,3,4),c(1,2,4,3),c(2,1,4,3),c(2,4,1,3),c(2,1,3,4),c(2,4,3,1),c(2,1,4,3),c(2,1,4,3),c(2,1,4,3),c(4,2,3,1)),nrow=10,byrow=T)
ratings <- matrix(c(c(0,0,2,3),c(0,0,2,4),c(0,1,2,2),c(0,2,3,3),c(0,1,2,3),c(0,1,2,4),c(0,2,4,3),c(1,1,2,3),c(1,0,3,4),c(1,1,3,4)),nrow=10,byrow=T)
M <- 4
sim2_alpha <- 0.15

# Visualize the ranking and rating data 
sim2_rankings_plot <- plot_rankings(rankings)
sim2_ratings_plot <- plot_ratings(ratings)

# Calculate original MB and new WMB maximum likelihood estimators
sim2_MLE_mb <- show_mb_estimate(rankings, ratings, M)
sim2_MLE_wmb <- show_wmb_estimate(rankings, ratings, M, alpha=sim2_alpha)

# Create plots of confidence intervals on p
sim2_CI_wmb <- ci_wmb(rankings=rankings,ratings=ratings,M=M,alpha=sim2_alpha,
                      interval=0.95,nsamples=200)
sim2_ci_plot <- plot_p_wmb(sim2_MLE_wmb, sim2_CI_wmb)

# Create plots showing how estimation results change with alpha
sim2_p_plot <- ASTAR_plot(rankings,ratings,M,n=101,parameter="p",return_plot=T)
sim2_theta_plot <- ASTAR_plot(rankings,ratings,M,n=101,parameter="theta",return_plot=T)

# Save all plots to files
png("conor/thesis/semester2/outputs/sim2_rankings_plot.png", width=6, height=5, units="in", res=1200)
sim2_rankings_plot
dev.off()
png("conor/thesis/semester2/outputs/sim2_ratings_plot.png", width=6, height=5, units="in", res=1200)
sim2_ratings_plot
dev.off()
png("conor/thesis/semester2/outputs/sim2_p_plot.png", width=6, height=5, units="in", res=1200)
sim2_p_plot
dev.off()
png("conor/thesis/semester2/outputs/sim2_theta_plot.png", width=6, height=5, units="in", res=1200)
sim2_theta_plot
dev.off()
png("conor/thesis/semester2/outputs/sim2_ci_plot.png", width=6, height=5, units="in", res=1200)
sim2_ci_plot
dev.off()


### TOY DATA 3 FROM RANKRATE PACKAGE

data("ToyData3")

rankings <- ToyData3$rankings
ratings <- ToyData3$ratings
M <- 4
toy3_alpha <- 0.15

# Visualize the ranking and rating data 
toy3_rankings_plot <- plot_rankings(rankings)
toy3_ratings_plot <- plot_ratings(ratings)

# Calculate original MB and new WMB maximum likelihood estimators
toy3_MLE_mb <- show_mb_estimate(rankings, ratings, M)
toy3_MLE_wmb <- show_wmb_estimate(rankings, ratings, M, alpha=toy3_alpha)

# Create plots of confidence intervals on p
toy3_CI_wmb <- ci_wmb(rankings=rankings,ratings=ratings,M=M,alpha=toy3_alpha,
                      interval=0.95,nsamples=200)
toy3_ci_plot <- plot_p_wmb(toy3_MLE_wmb, toy3_CI_wmb)

# Create plots showing how estimation results change with alpha
toy3_p_plot <- ASTAR_plot(rankings,ratings,M,n=101,parameter="p",return_plot=T)
toy3_theta_plot <- ASTAR_plot(rankings,ratings,M,n=101,parameter="theta",return_plot=T)

# Save all plots to files
png("conor/thesis/semester2/outputs/toy3_rankings_plot.png", width=6, height=5, units="in", res=1200)
toy3_rankings_plot
dev.off()
png("conor/thesis/semester2/outputs/toy3_ratings_plot.png", width=6, height=5, units="in", res=1200)
toy3_ratings_plot
dev.off()
png("conor/thesis/semester2/outputs/toy3_p_plot.png", width=6, height=5, units="in", res=1200)
toy3_p_plot
dev.off()
png("conor/thesis/semester2/outputs/toy3_theta_plot.png", width=6, height=5, units="in", res=1200)
toy3_theta_plot
dev.off()
png("conor/thesis/semester2/outputs/toy3_ci_plot.png", width=6, height=5, units="in", res=1200)
toy3_ci_plot
dev.off()


### SALMON DATA

M <- 4
salmon_alpha <- 0.15
# Monday ratings
likert_scale <- matrix(c(c(5,5,3,2),c(3,2,5,5),c(4,5,4,5),c(1,3,4,4),c(5,4,2,3),c(5,5,1,4),c(4,5,3,4),c(5,4,2,4),c(3,4,4,1),c(4,4,2,3)),nrow=10,byrow=T)
# Tuesday rankings
likert_scale2 <- matrix(c(c(4,4,3,2),c(5,4,3,4),c(5,4,3,3),c(3,3,1,2),c(5,4,3,4),c(5,5,5,4),c(4,4,3,4),c(3,4,2,4),c(5,5,5,5)),nrow=9,byrow=T)
ratings <- -1*((likert_scale-1)-4)
ratings2 <- -1*((likert_scale2-1)-4)
# Monday rankings
rankings <- matrix(c(c(3,4,1,2),c(3,4,1,2),c(3,1,4,2),c(3,1,4,2),c(1,4,2,3),c(2,1,4,3),c(3,4,1,2),c(3,1,4,2),c(3,1,2,4),c(3,1,4,2),c(2,4,1,3),c(4,3,1,2),c(3,2,1,4),c(3,4,1,2),c(3,4,1,2),c(2,1,4,3),c(3,4,1,2),c(4,3,1,2),c(3,4,2,1),c(3,2,4,1),c(3,1,4,2),c(3,1,2,4),c(4,3,1,2),c(3,4,2,1),c(3,4,2,1),c(3,4,1,2),c(3,4,2,1)),nrow=27,byrow=T)
# Tuesday rankings
rankings2 <- matrix(c(c(2,4,1,3),c(1,3,4,2),c(4,1,3,2),c(4,2,3,1),c(3,4,1,2),c(4,3,2,1),c(3,4,2,1),c(1,3,4,2),c(3,4,2,1),c(3,2,4,1),c(4,1,3,2),c(3,1,2,4),c(4,3,1,2),c(3,1,4,2),c(3,4,1,2),c(3,4,2,1),c(3,1,4,2),c(4,3,2,1),c(3,4,1,2),c(1,3,2,4),c(4,1,3,2),c(4,2,3,1),c(3,2,4,1),c(4,1,2,3),c(4,1,2,3),c(3,4,2,1),c(3,4,1,2)),nrow=27,byrow=T)
# Trim the rankings down to the same size as the ratings
# Selected seeds
# set.seed(325)
#set.seed(826678)
#set.seed(151627)
set.seed(297521)

# Seeds that are less accurate or interest versions of the above seeds
#set.seed(902647)
#set.seed(315204)
#set.seed(793020)
#set.seed(237737)
#set.seed(935252)
#set.seed(151627)
#set.seed(1498)

# Randomly sample the same number of rankings as there are available ratings
rankings <- rankings[sample(nrow(rankings),size=10,replace=FALSE),]
rankings2 <- rankings2[sample(nrow(rankings2),size=9,replace=FALSE),]

# Merge the two days, since the samples are from the same fillets
rankings <- rbind(rankings, rankings2)
ratings <- rbind(ratings, ratings2)

# Optional: Resample a smaller selection from the rankings and ratings
rankings <- rankings[sample(nrow(rankings),size=15,replace=TRUE),]
ratings <- ratings[sample(nrow(ratings),size=15,replace=TRUE),]

# Visualize the ranking and rating data 
salmon_rankings_plot <- plot_rankings(rankings)
salmon_ratings_plot <- plot_ratings(ratings)
  
# Calculate original MB and new WMB maximum likelihood estimators
salmon_MLE_mb <- show_mb_estimate(rankings, ratings, M)
salmon_MLE_wmb <- show_wmb_estimate(rankings, ratings, M, alpha=salmon_alpha)

# Create plots of confidence intervals on p
salmon_CI_wmb <- ci_wmb(rankings=rankings,ratings=ratings,M=M,alpha=salmon_alpha,
                      interval=0.95,nsamples=200)
salmon_ci_plot <- plot_p_wmb(salmon_MLE_wmb, salmon_CI_wmb)

# Create plots showing how estimation results change with alpha
salmon_p_plot <- ASTAR_plot(rankings,ratings,M,n=200,parameter="p",return_plot=T)
salmon_theta_plot <- ASTAR_plot(rankings,ratings,M,n=200,parameter="theta",return_plot=T)

# Save all plots to files
png("conor/thesis/semester2/outputs/salmon_rankings_plot.png", width=6, height=5, units="in", res=1200)
salmon_rankings_plot
dev.off()
png("conor/thesis/semester2/outputs/salmon_ratings_plot.png", width=6, height=5, units="in", res=1200)
salmon_ratings_plot
dev.off()
png("conor/thesis/semester2/outputs/salmon_p_plot.png", width=6, height=5, units="in", res=1200)
salmon_p_plot
dev.off()
png("conor/thesis/semester2/outputs/salmon_theta_plot.png", width=6, height=5, units="in", res=1200)
salmon_theta_plot
dev.off()
png("conor/thesis/semester2/outputs/salmon_ci_plot.png", width=6, height=5, units="in", res=1200)
salmon_ci_plot
dev.off()














