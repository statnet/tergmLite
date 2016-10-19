#--------------------------------------------------------------
#Load default parameters

primary_parameters  <- input_parameters_primary()
cd4_data            <- input_parameters_cd4_data()

#--- combine individual parameters into single list
evoparams <- c(primary_parameters, cd4_data)

#--------------------------------------------------------------
evoparams$nsims            = 1
evoparams$ncores=1

evoparams$popsumm_frequency = 25
evoparams$n_steps           = 365*50
evoparams$initial_pop       = 5000
evoparams$initial_infected  = .10*evoparams$initial_pop
evoparams$birth_model       = "poisson_birth_numbers"
evoparams$poisson_birth_lambda = 0.0137*(evoparams$initial_pop/100)
evoparams$trans_RR_age         = 1.0
evoparams$relation_dur  			 = 200
evoparams$trans_RR_insertive_anal_msm = 2.9  
evoparams$trans_RR_receptive_anal_msm = 17.3
evoparams$transmission_model  	      = "hughes"
evoparams$nw_form_terms = ~edges + concurrent + nodefactor("att1") + nodematch("att1", diff=F) + offset(nodematch('role', diff=TRUE, keep=1:2))  
evoparams$target_stats  <- c(0.7/2*evoparams$initial_pop, .10*evoparams$initial_pop, .80*(0.7*evoparams$initial_pop), .95*(0.7*evoparams$initial_pop/2))
evoparams$age_dist      <- seq(50, 10, -10/9)/1110
evoparams$generic_nodal_att_no_categories <- 2 
evoparams$generic_nodal_att_values <- 1:2 
evoparams$generic_nodal_att_values_props <- c(0.9, 0.1) 
evoparams$generic_nodal_att_values_props_births <- c(0.9, 0.1)

#add parameters that are functions of other input parameters
evoparams  <- input_parameters_derived(evoparams)

#convert raw parameter list into EpiModel object
evoparams <- do.call(EpiModel::param.net,evoparams)

#params to loop through----------------------------
#target_stats_vec <- "Run1"
#output name
#suffix=(rep("",1))
model_names= "Run11"

#estiamted_nw object names
suffix=(rep("",1))
nw_names = paste("evo_nw_",model_names,".RDATA",suffix,sep="")

#make sure param vecs same length

  load(file.path(outpath,nw_names))
  
  model_name = model_names
  #evoparams$target_stats   <- c(1750, target_stats_vec[ii])
  
  #--------------------------------------------------------------
  
  #-- create initial vector of infection status as an epimodel object
  infected_list <- EpiModel::init.net(i.num=evoparams$initial_infected,
                                      status.rand = FALSE)
  
  #--------------------------------------------------------------
  
  #---  Create list of modules to run for input into epimodel_control_fxn() below
  
  # ***   Note: initialize fxn must always be first and verbose fxn last AND death fxn
  # ***   must precede birth fxn (these conditions may change in future)
  # ***   treatment_fxn must be before update_vl and update_cd4
  
  evo_module_list<- list(
    "initialize.FUN"     = initialize_module,
    #    "plot_nw.FUN"        = plot_network_fxn,  
    "aging.FUN"          = vital_aging_module,
    "testing.FUN"        = social_testing_diagnosis_module,
    "treatment.FUN"      = social_treatment_module_john,
    "update_vl.FUN"      = viral_update_gamma,
    "update_cd4.FUN"     = viral_update_cd4_daily,
    "coital_acts.FUN"    = social_coital_acts_module,
    "trans.FUN"          = transmission_main_module,
    "trans_book.FUN"     = transmission_bookkeeping_module,
    "trans_cd4.FUN"      = transmission_cd4_module,
    "deaths.FUN"         = vital_deaths_module,
    "births.FUN"         = vital_births_module,
    "social_trans.FUN"   = social_attribute_transition_module,
    "summary.FUN"        = summary_module,
    "resim_nets.FUN"     = EpiModel::resim_nets,
    "verbose.FUN"        = NULL)
  
  
  #--- call epimodel's control fxn (load evonet modules into epimodel)
  evocontrol <- setup_epimodel_control_object(evonet_params = evoparams,
                                              module_list   = evo_module_list)
  
  #--------------------------------------------------------------
  
  if(isTRUE(hyak_par)){
    evomodel  <- EpiModelHPC::netsim_par(x = estimated_nw, 
                                         param = evoparams, 
                                         init = infected_list, 
                                         control = evocontrol)
  }else{
    evomodel  <- EpiModel::netsim(x = estimated_nw,
                                  param = evoparams,
                                  init = infected_list,
                                  control = evocontrol)
  }
  
  evomodel$epi <- NULL
  evomodel$stats <- NULL
  evomodel$control <- NULL
  
  
  plots_popsumm(evomodel,outpath=outpath,
                name=model_name,nw_stats=TRUE,max_points_rep=100,
                evoparams$popsumm_frequency)

  assign(model_name,evomodel)
  file_name <- paste(model_name,".RData",sep="")
  save(list=model_name,
       file = file.path(outpath,file_name) )
  remove(evomodel)
  
#end of loop

#--------------------------------------------------------------
