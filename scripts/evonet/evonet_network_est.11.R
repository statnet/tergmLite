#script to 1) estimate and 2) save to file a "network" object to use in evonet
#this step can not be done in hyak and must be done on linux and then object transferred
# to hyak
#
# NOTE! Very important!!!!
# "network" related parameters used to estimate network here need to be same on hyak run
# These "network" parameters include (if different from default)
# initital population size, target stats, relationship duration, role ,
# generic attribute, nw_form_terms, nw_coef_form, nw_constraints

#--------------------------------------------------------------
library(evonet)

#--------------------------------------------------------------
# change warnings to errors and open browser on error
options(warn=2) # turn warnings into error
options(error=browser) # go into debug mode on error

# change error setting back to default if desired
#options(error=NULL)
#options(warn=1)

#--------------------------------------------------------------
#Load default parameters

primary_parameters  <- input_parameters_primary()
cd4_data            <- input_parameters_cd4_data()

#--- combine individual parameters into single list
evoparams <- c(primary_parameters, cd4_data)
#--------------------------------------------------------------

evoparams$popsumm_frequency = 25
evoparams$n_steps           = 365*50
evoparams$initial_pop       = 5000
evoparams$initial_infected  = .10*evoparams$initial_pop
evoparams$birth_model       = "poisson_birth_numbers"
evoparams$poisson_birth_lambda = 0.0137*(evoparams$initial_pop/100)
evoparams$trans_RR_age         = 1.0
evoparams$relation_dur  			      	= 200
evoparams$trans_RR_insertive_anal_msm = 2.9  
evoparams$trans_RR_receptive_anal_msm = 17.3
evoparams$transmission_model  	      = "hughes"
evoparams$nw_form_terms = ~edges + concurrent + nodefactor("att1") + nodematch("att1", diff=F) + offset(nodematch('role', diff=TRUE, keep=1:2))  
evoparams$target_stats  <- c(0.7/2*evoparams$initial_pop, .10*evoparams$initial_pop, .80*(0.7*evoparams$initial_pop), .95*(0.7*evoparams$initial_pop/2))
evoparams$generic_nodal_att_no_categories <- 2 
evoparams$generic_nodal_att_values <- 1:2 
evoparams$generic_nodal_att_values_props <- c(0.9, 0.1) 
evoparams$generic_nodal_att_values_props_births <- c(0.9, 0.1)

modname = "Run11"

#add parameters that are functions of other input parameters
evoparams  <- input_parameters_derived(evoparams)

#convert raw parameter list into EpiModel object
evoparams <- do.call(EpiModel::param.net,evoparams)


#--------------------------------------------------------------
# initialize network

nw <- setup_initialize_network(evoparams)

#--------------------------------------------------------------

#run qaqc on input parameters
input_parameters_qaqc(evoparams)

#--------------------------------------------------------------

#estimate initial network (create argument list, then call fxn)
netest_arg_list <- list(
  nw            =  nw,
  formation     =  as.formula(evoparams$nw_form_terms),
  target.stats  =  evoparams$target_stats,
  coef.form     =  evoparams$nw_coef_form,
  constraints   =  as.formula(evoparams$nw_constraints),
  verbose       =  FALSE,
  coef.diss     =  dissolution_coefs( dissolution =  ~offset(edges),
                                      duration    =  evoparams$relation_dur,
                                      d.rate      =  3e-05) )

estimated_nw <- do.call(EpiModel::netest, netest_arg_list)

#--------------------------------------------------------------

save(estimated_nw,
     file = paste("evo_nw_",modname,".RDATA",sep=""))
remove(estimated_nw)

# change error setting back to default
options(error=NULL)
options(warn=1)
