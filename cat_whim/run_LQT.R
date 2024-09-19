#devtools::install_github('jdwor/LQT')
# Load the needed library
library(LQT)
library(doParallel)
library(stringr)

# Setup parallel
cl <- makeCluster(3)
registerDoParallel(cl)

########## Multi patients script ##########
# This script was run on Windows! Remember to change the BASE_DIR to your current project folder
########## Set up config structure ###########
BASE_DIR = 'C:/Users/leo_r/Desktop/cat_whim'
LQT_DIR = paste0(BASE_DIR, '/LQT')
MSK_DIR = paste0(LQT_DIR, '/lesion_masks')
RES_DIR = paste0(LQT_DIR, '/results')
LQT_DF_DIR = paste0(RES_DIR, '/dataframes')

if (! dir.exists(RES_DIR)){
  dir.create(RES_DIR)
}

if (! dir.exists(LQT_DIR)){
  dir.create(LQT_DIR)
}


lesion_paths = list.files(MSK_DIR,
                          full.names = TRUE)

pat_ids = str_extract(lesion_paths, "(?<=sub-)(.*)(?=_space)")

parcel_path = system.file("extdata","Other_Atlases",
                          "tpl-LQT_res-01_atlas-desikan_desc-dseg.nii.gz",package="LQT")

cfg = create_cfg_object(pat_ids=pat_ids,
                        lesion_paths=lesion_paths,
                        parcel_path=parcel_path,
                        out_path=RES_DIR,
)

########### Create Damage and Disconnection Measures ###########
start_time <- Sys.time()
# Get parcel damage for patients
get_parcel_damage(cfg, cores=1)
# Get tract SDC for patients
get_tract_discon(cfg, cores=1)
# Get parcel SDC and SSPL measures for patients
get_parcel_cons(cfg, cores=1)
end_time <- Sys.time()
end_time - start_time


data = compile_data(cfg, cores = 1)
list2env(data, .GlobalEnv);

if (! dir.exists(LQT_DF_DIR)){
  dir.create(LQT_DF_DIR)
}

######### Save Analysis-ready Datasets #############
write.csv(net.discon, paste0(LQT_DF_DIR, "/net_discon.csv"))
write.csv(net2net.discon, paste0(LQT_DF_DIR, "/net2net_discon.csv"))
write.csv(parc.damage, paste0(LQT_DF_DIR, "/parc_damage.csv"))
write.csv(parc.discon, paste0(LQT_DF_DIR, "/parc_discon.csv"))
write.csv(parc2parc.discon, paste0(LQT_DF_DIR, "/parc2parc_discon.csv"))
write.csv(tract.discon, paste0(LQT_DF_DIR, "/tract_discon.csv"))

######## Load and save SC matrix as csv ###########

for (subj in pat_ids){
  subj_dir = paste0(LQT_DIR, "/", subj)
  load(paste0(subj_dir, "/Parcel_Disconnection/", subj, "_tpl-LQT_res-01_atlas-desikan_desc-dseg_percent_parcel_mats.RData"))
  write.csv(pct_spared_sc_matrix, paste0(LQT_DIR, "/", subj, "/", "sub-", subj, "_pct_spared_sc_matrix.csv"))
  write.csv(pct_sdc_matrix, paste0(LQT_DIR, "/", subj, "/", "sub-", subj, "_pct_sdc_matrix.csv"))
}
