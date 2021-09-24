#projectname <- "ewing_ctDNA"
#samplename <- "409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
result_dir <- paste0(project_dir, "results_200821/")
in_path <- paste0(result_dir, "VAF_calculation/", samplename, "/Rdata/")

load(file.path(in_path, "VAFs_calculated.Rdata"))

result_dir <- paste0(project_dir, "results_200821/")
out_path <- paste0(result_dir, "VAF_calculation/", samplename, "/")