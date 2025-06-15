# set_model_parameters_v04.r

# v02 - change in normalization
# the change vs set_model_parameters_v02.r is in Processed/Total_counts/cells_born_counts.csv NOT in the parameters here.
# change v04 -> v05 in markov_premodel_prints_eps.r
#    pb.week <- c( 0, 1, 2, 4, 8, 12 )
#    if ( m$celltype == "B" ) {
#      pb.week <- c( 0, 1, 2, 4, 8 )
#    }
# the change v06 vs v05 is in Processed/Total_counts/cells_born_counts.csv -> NK NOT in the parameters here.
# the change v07 vs v06 is length of epsilon: 1hour -> 5 minutes
# the change v08 vs v07 is length of epsilon: 5 mins -> 2 mins

WEEK0 <- 0; 

HDR <- 0; 
U_HDR2 <- 1; 
U_WGHS <- 1; 

Q12_VALUE <- -0.5 
# birth rate FACTOR - multiple of total birth rates (e.g. 100000 for Treg) in .csv
BRF <- 1.1        
MAX_FLOW_TO_N3_FACTOR <- 10 

# # since v07 - 5 minutes
# EPS <- 9.920635e-05 * 5

# since v08 - 2 minutes
EPS <- 9.920635e-05 * 2


