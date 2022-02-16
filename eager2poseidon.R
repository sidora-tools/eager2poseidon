#!/usr/bin/env Rscript
if (!require('sidora.core')) {
  if(!require('remotes')) install.packages('remotes')
remotes::install_github('sidora-tools/sidora.core', quiet=T)
} else {library(sidora.core)}
if (!require(tidyverse)) {install.packages('tidyverse')}
require('dplyr', warn.conflicts=F)
require('readr')
require('purrr')
require('optparse')
if (!require('poseidonR')) {
  remotes::install_github('poseidon-framework/poseidonR', quiet=T)
}


## Validate genotype source option input
validate_geno_source <- function(option, opt_str, value, parser) {
  valid_entries=c("single", "double")
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\nInvalid genotype source: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

## Validate genotype ploidy option input
validate_ploidy <- function(option, opt_str, value, parser) {
  valid_entries=c("haploid", "diploid")
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\nInvalid genotype ploidy: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}


## Infer genetic sex from x-/y-rates. Caveat, no error bar taken into account.
infer_genetic_sex <- function(x_rate, y_rate) {
  case_when(
    `x_rate` > 0.70 & `y_rate` < 0.10 ~ "F",
    `x_rate` < 0.6 & `x_rate` > 0.3 & `y_rate` > 0.3 & `y_rate` < 0.6 ~ "M",
    TRUE ~ "U"
  )
}

## Calculate mean 5p1 damage for all shotgun libraries for which TF data exists, else all SG libs.
calculate_mean_5p1_damage <- function(ind, eager_data, lib_dmg) {
  # print(paste("ind:",ind))
  dmg <- NA
  ## Get individual's data
  x <- eager_data %>% 
    filter(grepl(ind, Ind))
  # print("Starting x")
  # print (x)
  
  ## If no SG data exists, return NA
  if (nrow(x %>% filter(grepl("SG",Seq))) == 0){
    # print(paste("Ind",ind,"had no SG data."))
    dmg <- NA
    warning(call.=F, 
            paste0("No SG data for individual ",ind,". DNA damage for individual in the janno file set to `n/a`."))
    ## If TF data exists, Only keep libraries with TF data that also have SG data.
  } else if (nrow(x %>% filter(grepl("TF",Seq))) > 0) {
    ## Create a tibble of presence absence for SG/TF data per library and seq_run.
    data_check <- x %>% select(Ind, Lib, Seq) %>% 
      filter(!is.na(Seq)) %>% 
      mutate(Seq=substr(Seq,1,2)) %>% 
      pivot_wider(names_from="Seq", values_from=Seq, values_fn=function(x) {! is.na(x)})
    # print(paste("Data check",ind))
    # print(data_check)
    ## Throw warning if TF data without SG exists
    if (nrow(data_check %>% filter( is.na(SG) & TF == TRUE)) >0) {
      bad_lib <- data_check %>% filter( is.na(SG) & TF == TRUE) %>% select(Ind,Lib)
      Lib_name=paste0(c(bad_lib,"TF"), collapse='.')
      warning(call.=F, 
              paste0("Detected TF data without matching SG for the same library. Data for '",
                     Lib_name,"' will not be included in 5p1 damage inference.\n",
                     "This could cause 'n/a's to be left in the final janno."))
      ##  Remove bad library from calculation
      x <- anti_join(x,bad_lib, by=c("Ind","Lib"))
    }
    ## Remove SG libraries that do not have TF data.
    x <- anti_join(x, data_check %>% filter(is.na(TF)) %>% select(Ind,Lib), by=c("Ind","Lib"))
    # print("Final x")
    # print (x)
    dmg <- x %>% filter(grepl("SG", Seq)) %>% inner_join(., lib_dmg, by=c("Ind", "Lib")) %>%
      summarise(total=sum(dedupped_reads.y), full_damage_5p1=damage_5p1.y * ( dedupped_reads.y / total)) %>%
      pull(full_damage_5p1) %>% sum()
    ## Return weighted mean dmg from all SG data
  } else {
    # print("Final x")
    # print (x)
    dmg <- x %>% filter(grepl("SG", Seq)) %>% inner_join(., lib_dmg, by=c("Ind", "Lib")) %>%
      summarise(total=sum(dedupped_reads.y), full_damage_5p1=damage_5p1.y * ( dedupped_reads.y / total)) %>%
      pull(full_damage_5p1) %>% sum()
  }
  return(dmg)
}

## Returns the max endogenous of all SG of all libraries for an individual, or NA if none are present
get_endogenous_data <- function(eager_data) {
  ## Only keep SG data
  x <- eager_data %>% 
    select(Ind,Lib,Seq,endogenous_dna) %>%
    filter(!is.na(endogenous_dna) & grepl("SG", Seq)) %>% 
    group_by(Ind) %>%
    summarise(Endogenous=max(endogenous_dna))
  return(x)
}

## Infer the correct calibration curve for quickcalibrate, based on the Pandora Latitude
select_calibration_curve <- function(Latitude) {
  case_when(
    Latitude >= 0 ~ "intcal20",
    Latitude < 0 ~ "shcal20",
    TRUE ~ "none"
  )
}

## Collect contamination data. All valid contamination estimates for each Library. TF when TF exists, else SG.
get_contamination_data <- function(eager_data) {
  x <- eager_data %>%
    filter(!is.na(x_contamination)) %>%
    filter(grepl("^TF", Seq)) %>%
    select(Ind, Lib, Seq, x_contamination, x_contamination_error)
  y <- eager_data %>%
    filter(!is.na(x_contamination)) %>%
    filter(grepl("^SG", Seq)) %>%
    anti_join(.,x, by=c("Ind", "Lib")) %>%
    select(Ind, Lib, Seq, x_contamination, x_contamination_error)
  
  contam <- bind_rows(x,y)
  multiple_libs <- contam %>% group_by(Ind) %>% summarise(n=n()) %>% filter(n>1)
  if (nrow(multiple_libs) > 1 ) {
    warning(call.=F, 
            paste0("More than one library with valid contamination data found for the following individuals:\n",
                   multiple_libs$Ind,"\nContamination of `n/a` will be added to janno for these individuals, as the field in the janno file\nshould reflect the estimate after merging all libraries of an individual."))
    contam <- anti_join(contam, multiple_libs, by=c("Ind"))
  }
  contam <- rename(contam, Xcontam=x_contamination, Xcontam_stderr=x_contamination_error)
  return(contam)
}

## Collect Pandora data at different levels
#### Site
get_site_info <- function(con, sites) {
  sites_table <- sidora.core::get_df(con, tab = "TAB_Site")
  inner_join(sites_table, sites, by=c("site.Site_Id"="Pandora_Site_Id")) %>%
    select(Group_Name=site.Full_Site_Id,
           Country=site.Country,
           Location=site.Locality,
           Site=site.Name,
           Latitude=site.Latitude,
           Longitude=site.Longitude)
}

#### Individual
get_individual_info <- function(con, individuals) {
  individual_table <- sidora.core::get_df(con, tab = "TAB_Individual")
  inner_join(individual_table, individuals, by=c("individual.Full_Individual_Id"="Ind")) %>%
    select(Individual_ID=individual.Full_Individual_Id, 
           Collection_ID=individual.Archaeological_ID)
}

#### Sample
get_sample_info <- function(con, samples) {
  sample_table <- sidora.core::get_df(con, tab = "TAB_Sample") %>% convert_all_ids_to_values(., con = con)
  inner_join(sample_table, samples, by=c("sample.Full_Sample_Id"="Pandora_Sample_Id")) %>%
    select(Sample_Id=sample.Full_Sample_Id,
           sample.Type_Group, 
           sample.Type) %>%
    unite(col="Source_Tissue", sep=":", na.rm=T, tolower(sample.Type_Group), tolower(sample.Type))
}

#### Library
get_library_info <- function(con, libraries) {
  library_table <- sidora.core::get_df(con, tab = "TAB_Library") %>% convert_all_ids_to_values(., con = con)
  inner_join(library_table, libraries, by=c("library.Full_Library_Id"="Pandora_Library_Id")) %>%
    select(Library_Id=library.Full_Library_Id, 
           Protocol=library.Protocol)
}

## Infer row type by ID length and TF/SG
# infer_data_type <- function()
  
## MAIN ##

## Parse arguments ----------------------------
parser <- OptionParser(usage = "%prog [options] .credentials")
parser <- add_option(parser, c("-i", "--input"), type = 'character', 
                     action = "store", dest = "input", 
                     help = "A path to the eager output 
                    directory to turn into a poseidon package.")
parser <- add_option(parser, c("-n", "--outPackageName"), type = 'character',
                     action = "store", dest = "package_name", 
                     help = "The output package name")
parser <- add_option(parser, c("-o", "--outPackagePath"), type = 'character',
                     action = "store", dest = "package_path", 
                     help = "The output package directory path")
parser <- add_option(parser, c("-g", "--genotypeSource"), type = 'character',
                     action = "callback", dest = "genotype_source", 
                     metavar="source",
                     callback = validate_geno_source,
                     help = "If both single and double stranded genotypes exist in the eager output, which set should be converted into a poseidon package?\n\t\tMerging the two (in case some individuals only exist in either of the two) is currently not supported.")
parser <- add_option(parser, c("-p", "--genotypePloidy"), type = 'character',
                     action = "callback", dest = "genotype_ploidy",
                     metavar="ploidy",
                     callback = validate_ploidy, default="haploid",
                     help = "The genotype ploidy of the genotypes produced by eager. This value will be used as the 'Genotype_Ploidy' of all individuals in\n\t\tthe output janno file. The default value for this parameter is 'haploid'.")
arguments <- parse_args(parser, positional_arguments = 1)

opts <- arguments$options
cred_file <- arguments$args

###########
## EAGER ##
###########

## Load eager output data
raw_eager_data <- read_tsv(paste0(opts$input, "/multiqc/multiqc_data/multiqc_general_stats.txt"), na = c("", "N/A", "NA")) 
eager_data <- raw_eager_data %>% 
  ## Only Keep columns of interest
  select(
    `Sample`, 
    rate_x=`SexDetErrmine_mqc-generalstats-sexdeterrmine-RateX`,
    rate_y=`SexDetErrmine_mqc-generalstats-sexdeterrmine-RateY`,
    endogenous_dna=`endorSpy_mqc-generalstats-endorspy-endogenous_dna`,
    convered_snps=`snp_coverage_mqc-generalstats-snp_coverage-Covered_Snps`,
    total_snps=`snp_coverage_mqc-generalstats-snp_coverage-Total_Snps`, 
    damage_5p1=`DamageProfiler_mqc-generalstats-damageprofiler-5_Prime1`,
    x_contamination_snps=`nuclear_contamination_mqc-generalstats-nuclear_contamination-Num_SNPs`,
    x_contamination=`nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_estimate`,
    x_contamination_error=`nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_SE`,
    duplication_rate=`Picard_mqc-generalstats-picard-PERCENT_DUPLICATION`,
    filtered_mapped_reads=`Samtools Flagstat (post-samtools filter)_mqc-generalstats-samtools_flagstat_post_samtools_filter-mapped_passed`
  ) %>% 
  ## Remove rows with all but the sample name being NA (sequencing lanes)
  filter(rowSums(is.na(.)) != ncol(.) - 1) %>%
  ## Calculate dedupped reads. Cannot use Qualimap because it is calculated after merging TF and SG data.
  mutate(dedupped_reads=round(filtered_mapped_reads * ( 1 - duplication_rate ),0)) %>%
  separate(Sample, into=c("Ind","Lib","Seq"), sep="\\.", fill='right') %>% 
  ##Set Xcont to 0 if negative, NA if over 1
  mutate(
    Pandora_Site_Id=substr(Ind,1,3),
    Pandora_Library_Id=paste(sep=".", Ind,Lib),
    Pandora_Sample_Id=substr(Pandora_Library_Id,1,8),
    x_contamination=case_when(
      x_contamination < 0 ~ 0,  ## If estimate is negative, replace with 0 as per janno format specs.
      x_contamination_snps != 0 ~ x_contamination
      ## When x_contamination_snps is 0, add NA
    ), Genetic_Sex=infer_genetic_sex(rate_x, rate_y)
  )

## Create tibble of SG damage_5p1 per Library
lib_dmg <- eager_data %>% select(Ind,Lib,Seq, damage_5p1, dedupped_reads) %>% 
  filter(grepl("SG", Seq))

## Get contamination data
contamination_data <- get_contamination_data(eager_data)

## Add damage_calculation to eager-janno
eager_janno <- eager_data %>% 
  filter(is.na(Lib)) %>% ## Keep only individual level entries
  mutate(damage_5p1=map_dbl(Ind, calculate_mean_5p1_damage, eager_data, lib_dmg)) %>%
  select(-starts_with("x_cont"), -Lib, -Seq, -endogenous_dna,-duplication_rate, -filtered_mapped_reads, -dedupped_reads) %>%
  ## Add nuclear contamination 
  left_join(., contamination_data %>% select(-Lib, -Seq), by=c("Ind")) %>%
  mutate(
    Xcontam=ifelse(Genetic_Sex == "F", NA, Xcontam),
    Xcontam_stderr=ifelse(Genetic_Sex == "F", NA, Xcontam_stderr)
  ) %>%
  left_join(., get_endogenous_data(eager_data), by="Ind") %>%
  ## Add genotype ploidy
  mutate(Genotype_Ploidy=opts$genotype_ploidy)
eager_janno
cat("Inference from eager results complete.")

#############
## PANDORA ##
#############

con <- get_pandora_connection(cred_file)

## Get complete pandora table
complete_pandora_table <- join_pandora_tables(
  get_df_list(
    c(make_complete_table_list(
      c("TAB_Site", "TAB_Raw_Data")
    )), con = con
  )
) %>% convert_all_ids_to_values(., con = con)


results <- inner_join(complete_pandora_table, eager_janno %>% select(Ind), by=c("individual.Full_Individual_Id"="Ind"))
  x %>% select(Individual_ID=individual.Full_Individual_Id, 
               Collection_ID=individual.Archaeological_ID, 
               Group_Name=site.Full_Site_Id,
               Country=site.Country,
               Location=site.Locality,
               Site=site.Name,
               Latitude=site.Latitude,
               Longitude=site.Longitude,
               Date_C14_Uncal_BP=individual.C14_Uncalibrated,
               Date_C14_Uncal_BP_Err=individual.C14_Uncalibrated_Variation,
               Date_C14_Labnr=individual.C14_Id,
               Dating_Note=individual.C14_Info, ## Keep note to make issues easier to spot.
               Calibration_Curve=select_calibration_curve(site.Longitude)
               ) %>%
    mutate(Location=if_else(Location == "", "n/a", Location),
           Date_BC_AD_Median = quickcalibrate(Date_C14_Uncal_BP, Date_C14_Uncal_BP_Err,Calibration_Curve, allowOutside=T) %>% .$Date_BC_AD_Median,
           Date_BC_AD_Start = quickcalibrate(Date_C14_Uncal_BP, Date_C14_Uncal_BP_Err, Calibration_Curve, allowOutside=T)%>% .$Date_BC_AD_Start,
           Date_BC_AD_Stop = quickcalibrate(Date_C14_Uncal_BP, Date_C14_Uncal_BP_Err, Calibration_Curve, allowOutside=T) %>% .$Date_BC_AD_Stop
    ) %>% View()

## Remove Sites that have no individuals yet
complete_pandora_table %>% filter(!is.na(individual.Full_Individual_Id)) %>%
  

