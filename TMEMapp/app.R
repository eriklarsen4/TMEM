
library(shiny)
library(bslib)
library(TMEM)
library(Endo)
library(forcats)


# import datasets----

"%notin%" <- Negate("%in%")

aDRG_DEA_results <- TMEM::aDRG_DEA_results
aDRG_TPM <- TMEM::aDRG_TPM

read.table('C:/Users/Erik/Desktop/Programming/R/Bio/Itch/aDRG_GO_analysis_BP.txt', )

e13DRG_DEA_results <- read.csv('C:/Users/Erik/Downloads/DESeq2_result_file_on_GT_e13_PARA.csv', header = F)
p0DRG_DEA_results <- read.csv('C:/Users/Erik/Downloads/DESeq2_result_file_on_GT_p0_mu234_wtBCD_PARA.csv', header = T)
p10DRG_DEA_results <- read.csv('C:/Users/Erik/Downloads/DESeq2_result_file_on_GT_p10_mu123_wt1CE_PARA.csv', header = T)

eDRG_TPM <- read.csv('C:/Users/Erik/Downloads/DESeq2_Normalized_Counts_Entire_Timecourse_by_Genotype.csv', header = F)

CRAPome_results <- Endo::CRAPome_results
IPMS_results <- Endo::IPMS_counts

Proteomics_GO_BP_results <- Endo::GO_BP_results
Proteomics_GO_CC_results <- Endo::GO_CC_results
Proteomics_GO_MF_results <- Endo::GO_MF_results

# RNA-seq results wrangling ----

aDRG_DEA_results <- aDRG_DEA_results |>
  dplyr::mutate(`%WT` = as.numeric(`%WT`)*100)

e13DRG_DEA_results <- e13DRG_DEA_results |>
  dplyr::mutate(`%WT` = 2^(V3)*100)

p0DRG_DEA_results <- p0DRG_DEA_results |>
  dplyr::mutate(`%WT` = 2^(log2.FC.)*100)

p10DRG_DEA_results <- p10DRG_DEA_results |>
  dplyr::mutate(`%WT` = 2^(log2.FC.)*100)

colnames(e13DRG_DEA_results) = colnames(aDRG_DEA_results)
colnames(p0DRG_DEA_results) = colnames(aDRG_DEA_results)
colnames(p10DRG_DEA_results) = colnames(aDRG_DEA_results)


itch_GO_info <- TMEM::get_GO_info(list_of_interest = aDRG_DEA_results |> dplyr::filter(AdjP <= 0.05) |> dplyr::distinct(GeneID) |> unlist() |> as.character(), species = 'mouse')

aDRG_DEA_results |>
  dplyr::inner_join(itch_GO_info$GO_info_by_term_df)

itch_GO_info$unique_GOs

itch_GO_info$GO_info_by_term_df

itch_GO_info |>
  dplyr::mutate(fold.Enrichment = as.numeric(fold.Enrichment)) |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |>
  dplyr::filter(!grepl(GO_Term, pattern = 'Unclassified', ignore.case = T)) |>
  dplyr::mutate(sig = dplyr::case_when(

    FDR <= 0.001 ~ "***",
    dplyr::between(FDR, 0.001, 0.0499) ~ "**",
    dplyr::between(FDR, 0.01, 0.05) ~ "*",
    TRUE ~ "NS"), .after = FDR
  )

GO_info_by_term_df_sig <- Proteomics_GO_results |>
  dplyr::mutate(db = 'geneontology.org') |>
  dplyr::full_join(x$GO_info_by_term_df |>
                     dplyr::mutate(db = 'bioconductor')) |>
  dplyr::filter(!is.na(GO_Term)) |>

  dplyr::mutate(fold.Enrichment = as.numeric(fold.Enrichment)) |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |>
  dplyr::filter(!grepl(GO_Term, pattern = "Unclassified", ignore.case = T)) |>
  dplyr::mutate(sig = dplyr::case_when(

    FDR <= 0.001 ~ "***",
    dplyr::between(FDR, 0.001, 0.0499) ~ "**",
    dplyr::between(FDR, 0.01, 0.05) ~ "*",
    TRUE ~ "NS"), .after = FDR
  )
rm(x)

# Proteomics results wrangling ----
IPMS_results <- IPMS_results |>
  dplyr::select(-1) |>
  dplyr::rename_with(

    ~IPMS_results |>
      dplyr::select(-1) |>
      dplyr::slice(1) |>
      unlist() |>
      as.character()
  ) |>
  dplyr::slice(-1) |>
  dplyr::slice_tail(n = -1) |>
  dplyr::select(-dplyr::contains("307-411")) |>
  dplyr::rename(
    Proteins.Identified = dplyr::contains("Proteins")
  ) |>
  dplyr::rename(
    Accession.Number = dplyr::contains("Accession")
  ) |>
  dplyr::rename(
    MW = dplyr::contains("Molecular")
  ) |>
  dplyr::filter(
    !grepl(Accession.Number, pattern = 'DECOY')
  ) |>
  dplyr::rename_with(
    ~gsub(.x, pattern = '[0-9]+_', replacement = '')
  ) |>

  dplyr::rename_with(
    ~gsub(.x, pattern = '^m', replacement = '')
  ) |>

  dplyr::rename(BirAV5_1 = `BirA-V1`,
                BirAV5_2 = `BirA-V2`,
                TMEMV5_1 = `TMEM-V1`,
                TMEMV5_2 = `TMEM-V2`
  ) |>
  dplyr::mutate(
    HumanGeneID = stringr::str_extract(
      Proteins.Identified,
      pattern = 'GN=[:alnum:]+') |>
      gsub(Proteins.Identified, pattern = 'GN=', replacement = ''),
    .after = Proteins.Identified
  ) |>

  dplyr::mutate(
    Proteins.Identified = gsub(Proteins.Identified,
                               pattern = ' O.+$',
                               replacement = '')
  ) |>

  dplyr::mutate(
    MW = gsub(MW, pattern = ' kDa', replacement = '') |> as.numeric()
  ) |>

  dplyr::mutate(
    dplyr::across(dplyr::contains("_"), ~dplyr::if_else(.x == '?', '0', .x)),
    dplyr::across(dplyr::contains("_"), ~dplyr::if_else(.x == '', '0', .x)),
    dplyr::across(dplyr::contains("_"), ~as.numeric(.x)
    ),
    dplyr::across(dplyr::contains("_"), ~dplyr::if_else(.x == 0, 0.1, .x)
    ),
    dplyr::across(dplyr::contains("_"), .fns = list("MW_norm" = ~.x/MW)
    )
  ) |>
  dplyr::mutate(

    V5ControlAvg = (BirAV5_1_MW_norm + BirAV5_2_MW_norm)/2,
    mycControlAvg = (GFPmyc_1_MW_norm + GFPmyc_2_MW_norm)/2,
    V5TMEMAvg = (TMEMV5_1_MW_norm + TMEMV5_2_MW_norm)/2,
    mycTMEMAvg = (TMEMmyc_1_MW_norm + TMEMmyc_2_MW_norm)/2,

    avgV5_FC = V5TMEMAvg/V5ControlAvg,
    avgmyc_FC = mycTMEMAvg/mycControlAvg,

    ControlAvg = (BirAV5_1_MW_norm+
                    BirAV5_2_MW_norm+
                    GFPmyc_1_MW_norm+
                    GFPmyc_2_MW_norm
    )/4,

    TMEMAvg = (TMEMV5_1_MW_norm+
                 TMEMV5_2_MW_norm+
                 TMEMmyc_1_MW_norm+
                 TMEMmyc_2_MW_norm
    )/4,

    avgFC = TMEMAvg/ControlAvg

  ) |>
  dplyr::mutate(

    MouseGeneID = paste(substr(HumanGeneID, 1, 1),
                        tolower(substr(HumanGeneID, 2, 10)), sep = ''),

    MouseGeneID = dplyr::case_when(
      grepl(Accession.Number, pattern = 'p53_HUMAN') ~ "Trp53",
      HumanGeneID == 'birA' ~ HumanGeneID,
      Proteins.Identified == 'Green fluorescent protein' ~ HumanGeneID,
      TRUE ~ MouseGeneID)
  ) |>

  dplyr::relocate(MouseGeneID, .after = HumanGeneID) |>
  dplyr::arrange(desc(avgFC))

Compiled_Candidate_List <- IPMS_results |>
  dplyr:::arrange(desc(avgFC)) |>
  dplyr::filter(
    avgFC > 2,
    avgV5_FC > 2,
    avgmyc_FC > 2,


    (TMEMV5_1 & TMEMV5_2 > 5),

    (TMEMmyc_1 & TMEMmyc_2 > 5),

    (
      TMEMV5_1 > 10 |
        TMEMV5_2 > 10 |
        TMEMmyc_1 > 10 |
        TMEMmyc_2 > 10
    )

  ) |>


  dplyr::filter(

    !grepl(
      HumanGeneID,
      pattern = 'AHNAK|ALB|ANXA(1|2)|ASS1|ATP5|CFL1|CLTC|CSDA|EDARRAD|ENO|FARS(A|B)|FASN|FKSG30|GAPDH|GFAP|ILF(2|3)|LDHA|LGALS(1|3)|PCB(1|2)|PDIA6|PHB|PHB2|POTE2|PPIA|PRDX(1|2|3|4)|SERPIN(H1|A11)|SLC25A|TUFM|TXN|UBA52')

  ) |>


  dplyr::filter(

    !grepl(
      HumanGeneID,
      pattern = 'EEF1|EIF|HNRN|HSPA|KRT|TUB(A|B)|RPL|HIST2|H2A|RRP|(X|I)PO|NUP|DNAJ')

  ) |>

  dplyr::filter(HumanGeneID != "TMEM184B") |>
  dplyr::select(HumanGeneID) |>
  unlist() |> as.character()


Crapome_results <- CRAPome_results |>
  dplyr::filter(!grepl(Mapped, pattern = '/|identifier')) |>
  dplyr::select(1:7) |>
  dplyr::rename(User_input = User,
                Mapped_Gene_Symbol = Input,
                Num_Exp_Found = Mapped,
                Num_Exp_Total = Symbol,
                Avg_SC = Num,
                Max_SC = of) |>
  dplyr::mutate(dplyr::across(c(3,5:7), ~as.numeric(.x))) |>
  dplyr::select(-Gene)
rm(CRAPome_results)


Crapome_hits <- Crapome_results |>
  dplyr::filter(Num_Exp_Found <= 35.8) |> # fewer than 5% of experiments
  dplyr::select(Mapped_Gene_Symbol) |>
  unlist() |>
  as.character()


IPMS_results_filtered <- IPMS_results |>
  dplyr::mutate(poi = 'Background', .after = MouseGeneID) |>
  dplyr::full_join(Crapome_hits |>
                     as.data.frame() |>
                     dplyr::mutate(HumanGeneID = Crapome_hits)
  ) |>
  dplyr::mutate(poi = dplyr::case_when(

    !is.na(Crapome_hits) ~ 'Candidates',
    HumanGeneID == 'TMEM184B' ~ 'TMEM184B',
    TRUE ~ poi)

  ) |>
  dplyr::select(-Crapome_hits)

  # Proteomics GO wrangling ----
Proteomics_GO_BP_results <- Proteomics_GO_BP_results |>
  dplyr::rename_with(
    ~gsub(.x, pattern = '\\.\\.', replacement = '')
  ) |>

  dplyr::rename_with(
    ~gsub(.x, pattern = 'Homo.sapiens.|upload_1|\\.$', replacement = '')
  ) |>

  dplyr::rename(
    GO_Term = GO.biological.process.complete
  ) |>

  dplyr::rename(
    GO_Term_Size = REFLIST20580
  ) |>

  dplyr::rename(
    Overlap = `138`
  ) |>
  dplyr::mutate(
    dplyr::across(c(2:4,7,8), ~as.numeric(.x))
  ) |>


  tidyr::separate(
    GO_Term, into = c("GO_Term", "GO_Term_ID"), sep = '\\s\\(GO:'
  ) |>

  dplyr::mutate(
    GO_Term_ID = paste0('GO:',
                        gsub(GO_Term_ID, pattern = '\\)', replacement = ''))
  ) |>

  dplyr::mutate(GO_Type = "CC",
                .after = GO_Term_ID
  ) |>

  dplyr::mutate(
    fold.Enrichment = stringr::str_extract(fold.Enrichment,
                                           pattern = '[:alnum:]+') |>
      as.numeric(fold.Enrichment)
  ) |>

  dplyr::arrange(
    FDR, desc(fold.Enrichment)
  )

Proteomics_GO_CC_results <- Proteomics_GO_CC_results |>
  dplyr::rename_with(
    ~gsub(.x, pattern = '\\.\\.', replacement = '')
  ) |>

  dplyr::rename_with(
    ~gsub(.x, pattern = 'Homo.sapiens.|upload_1|\\.$', replacement = '')
  ) |>

  dplyr::rename(
    GO_Term = GO.cellular.component.complete
  ) |>

  dplyr::rename(
    GO_Term_Size = REFLIST20580
  ) |>

  dplyr::rename(
    Overlap = `138`
  ) |>
  dplyr::mutate(
    dplyr::across(c(2:4,7,8), ~as.numeric(.x))
  ) |>

  tidyr::separate(
    GO_Term, into = c("GO_Term", "GO_Term_ID"), sep = '\\s\\(GO:'
  ) |>

  dplyr::mutate(
    GO_Term_ID = paste0('GO:',
                        gsub(GO_Term_ID, pattern = '\\)', replacement = ''))
  ) |>

  dplyr::mutate(GO_Type = "CC",
                .after = GO_Term_ID
  ) |>

  dplyr::mutate(
    fold.Enrichment = stringr::str_extract(fold.Enrichment,
                                           pattern = '[:alnum:]+') |>
      as.numeric(fold.Enrichment)
  ) |>

  dplyr::arrange(
    FDR, desc(fold.Enrichment)
  )

Proteomics_GO_MF_results <- Proteomics_GO_MF_results |>
  dplyr::rename_with(
    ~gsub(.x, pattern = '\\.\\.', replacement = '')
  ) |>

  dplyr::rename_with(
    ~gsub(.x, pattern = 'Homo.sapiens.|upload_1|\\.$', replacement = '')
  ) |>

  dplyr::rename(
    GO_Term = GO.molecular.function.complete
  ) |>

  dplyr::rename(
    GO_Term_Size = REFLIST20580
  ) |>

  dplyr::rename(
    Overlap = `138`
  ) |>

  dplyr::mutate(
    dplyr::across(c(2:4,7,8), ~as.numeric(.x))
  ) |>

  tidyr::separate(
    GO_Term, into = c("GO_Term", "GO_Term_ID"), sep = '\\s\\(GO:'
  ) |>

  dplyr::mutate(
    GO_Term_ID = paste0('GO:',
                        gsub(GO_Term_ID, pattern = '\\)', replacement = ''))
  ) |>
  dplyr::mutate(
    GO_Type = "MF", .after = GO_Term_ID
  ) |>

  dplyr::mutate(
    fold.Enrichment = stringr::str_extract(fold.Enrichment,
                                           pattern = '[:alnum:]+') |>
      as.numeric(fold.Enrichment)
  ) |>

  dplyr::arrange(
    FDR, desc(fold.Enrichment)
  )

Proteomics_GO_results <- Proteomics_GO_BP_results |>
  dplyr::full_join(Proteomics_GO_CC_results) |>
  dplyr::full_join(Proteomics_GO_MF_results) |>
  dplyr::filter(!grepl(GO_Term, pattern = 'unclassified', ignore.case = T))
rm(list = ls()[which(grepl(ls(), pattern = 'Proteomics_GO_(BP|MF|CC)_results')== T)])

x <- TMEM::get_GO_info(list_of_interest = Crapome_hits, species = 'human')

GO_info_by_term_df_sig <- Proteomics_GO_results |>
  dplyr::mutate(db = 'geneontology.org') |>
  dplyr::full_join(x$GO_info_by_term_df |>
                     dplyr::mutate(db = 'bioconductor')) |>
  dplyr::filter(!is.na(GO_Term)) |>

  dplyr::mutate(fold.Enrichment = as.numeric(fold.Enrichment)) |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |>
  dplyr::filter(!grepl(GO_Term, pattern = "Unclassified", ignore.case = T)) |>
  dplyr::mutate(sig = dplyr::case_when(

    FDR <= 0.001 ~ "***",
    dplyr::between(FDR, 0.001, 0.0499) ~ "**",
    dplyr::between(FDR, 0.01, 0.05) ~ "*",
    TRUE ~ "NS"), .after = FDR
  )
rm(x)

IPMS_results_filtered_full <- IPMS_results_filtered |>
  dplyr::mutate(

    poi = dplyr::case_when(
      grepl(Proteins.Identified, pattern = 'Green fluorescent') ~ 'Control myc',
      grepl(Proteins.Identified, pattern = 'BirA') ~ 'Control V5',
      TRUE ~ poi)

  )

# Define UI for application that draws a histogram ----
ui <- shiny::navbarPage(

  title = 'Bhattacharya Lab NGS datasets',

  theme = bslib::bs_theme(
    preset = 'cosmo',
    bg = 'white',
    fg = '#000080',
    primary = '#a40000',
    secondary = '#000080',
    sucesss = '#a40000'
  ),

  collapsible = TRUE,
  inverse = FALSE,

  # shiny::navbarMenu("NGS Dataset")

  shiny::tabPanel(strong("RNA-seq"),
                  shiny::navbarMenu(strong('Single Dataset'),
                                    shiny::tabPanel('Adult Dorsal Root Ganglia',
                                                    shiny::navbarMenu(
                                                      strong('Volcano'),
                                                      shiny::tabPanel(
                                                        strong('Volcano'),
                                                        shiny::sidebarLayout(
                                                          sidebarPanel = shiny::sidebarPanel(
                                                            # shiny::selectizeInput(inputId = 'VizChoice',
                                                            #                       label = p(strong('Visualization Type'), br('(Choose One)')),
                                                            #                       choices = c('Volcano', 'Heatmap', 'MA Plot', 'Differentially Expressed Genes Table', 'Normalized Counts Table')),

                                                            shiny::sliderInput(inputId = 'SignificanceSlider',
                                                                               label = 'Adjusted P-value',
                                                                               min = aDRG_DEA_results |> dplyr::slice_head() |> dplyr::select(AdjP) |> unlist() |> as.numeric(),
                                                                               max = aDRG_DEA_results |> dplyr::slice_tail() |> dplyr::select(AdjP) |> unlist() |> as.numeric(),
                                                                               value = c(aDRG_DEA_results |> dplyr::slice_head() |> dplyr::select(AdjP) |> unlist() |> as.numeric(),
                                                                                         aDRG_DEA_results |> dplyr::slice_tail() |> dplyr::select(AdjP) |> unlist() |> as.numeric()),
                                                                               dragRange = TRUE),

                                                            shiny::sliderInput(inputId = 'FCSlider',
                                                                               label = 'Expression FC',
                                                                               min = aDRG_DEA_results |> dplyr::arrange(log2.FC.) |> dplyr::slice_head() |> dplyr::select(log2.FC.) |> unlist() |> as.numeric(),
                                                                               max = aDRG_DEA_results |> dplyr::arrange(log2.FC.) |> dplyr::slice_tail() |> dplyr::select(log2.FC.) |> unlist() |> as.numeric(),
                                                                               value = c(aDRG_DEA_results |> dplyr::arrange(log2.FC.) |> dplyr::slice_head() |> dplyr::select(log2.FC.) |> unlist() |> as.numeric(),
                                                                                         aDRG_DEA_results |> dplyr::arrange(log2.FC.) |> dplyr::slice_tail() |> dplyr::select(log2.FC.) |> unlist() |> as.numeric()),
                                                                               dragRange = TRUE),

                                                            shiny::checkboxInput(inputId = 'FilterBySignificance',
                                                                                 label = p(strong('Significantly Differentially Expressed')),
                                                                                 value = FALSE),

                                                            shiny::selectizeInput(inputId = 'SearchByGene',
                                                                                  label = p(strong('Filter By Gene'), br('(Choose Any)')),
                                                                                  choices = c('All', c(aDRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character())),
                                                                                  selected = 'All',
                                                                                  multiple = TRUE),

                                                            shiny::selectizeInput(inputId = 'MapAnyGOTerm',
                                                                                  label = p(strong('Map Any GO Term'), br('(Choose One)')),
                                                                                  choices = c('None', itch_GO_info$unique_GOs)),

                                                            shiny::selectizeInput(inputId = 'MapEnrichedGOTerm',
                                                                                  label = p(strong('Map Significantly Enriched GO Term'), br('(Choose One)')),
                                                                                  choices = c('None', itch_GO_info$))
                                                            ),

                                                            mainPanel = shiny::mainPanel(
                                                              shiny::plotOutput('DRGVolcano', inline = T)
                                                            )
                                                        )
                                                      ),

                                                      strong('Heatmap'),
                                                      shiny::tabPanel(
                                                        strong('Heatmap'),
                                                        shiny::sidebarLayout(
                                                          sidebarPanel = shiny::sidebarPanel(


                                                            shiny::selectizeInput(inputId = 'SearchByGene',
                                                                                  label = p(strong('Filter By Gene'), br('(Choose Any)')),
                                                                                  choices = c('All', c(aDRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character())),
                                                                                  selected = 'All',
                                                                                  multiple = TRUE),

                                                            shiny::checkboxInput(inputId = 'FilterBySignificance',
                                                                                 label = p(strong('Significantly Differentially Expressed')),
                                                                                 value = FALSE)
                                                          ),

                                                          mainPanel = shiny::mainPanel(
                                                            shiny::plotOutput('DRGHeatmap', inline = T)
                                                          )
                                                        )
                                                      ),

                                                      strong('MA Plot'),
                                                      strong('Differentially Expressed Genes Table'),
                                                      strong('Normalized Counts Table'),

                                                    ),
                                    shiny::navbarMenu('Embryonic Dorsal Root Ganglia',
                                                      shiny::navbarMenu('Comparison',
                                                                        shiny::tabPanel('Age',
                                                                                        shiny::sidebarLayout(
                                                                                          sidebarPanel = shiny::sidebarPanel(
                                                                                            shiny::selectizeInput(inputId = 'VizChoice',
                                                                                                                  label = p(strong('Visualization Type'), br('(Choose One)')),
                                                                                                                  choices = c('Volcano', 'Heatmap', 'MA Plot', 'Differentially Expressed Genes Table', 'Normalized Counts Table')),
                                                                                            shiny::checkboxInput(inputId = 'FilterBySignificance',
                                                                                                                 label = p(strong('Significantly Differentially Expressed')),
                                                                                                                 value = FALSE),
                                                                                            shiny::selectizeInput(inputId = 'FilterByDirection',
                                                                                                                  label = p(strong('Regulation'), br('(Choose One)')),
                                                                                                                  choices = c('Either', 'Up', 'Down'),
                                                                                                                  selected = 'Either',
                                                                                                                  multiple = FALSE,
                                                                                                                  options = list(maxItems = 1)),
                                                                                            shiny::selectizeInput(inputId = 'SearchByGene',
                                                                                                                  label = p(strong('Filter By Gene'), br('(Choose Any)')),
                                                                                                                  choices = c('All',
                                                                                                                              c(e13DRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character(),
                                                                                                                                p0DRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character(),
                                                                                                                                p10DRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character()) |> unique() ),
                                                                                                                  selected = 'All',
                                                                                                                  multiple = TRUE)
                                                                                          )
                                                                                        )),
                                                                        shiny::tabPanel('Genotype',
                                                                                        shiny::sidebarLayout(
                                                                                          sidebarPanel = shiny::sidebarPanel(
                                                                                            shiny::selectizeInput(inputId = 'VizChoice',
                                                                                                                  label = p(strong('Visualization Type'), br('(Choose One)')),
                                                                                                                  choices = c('Volcano', 'Heatmap', 'MA Plot', 'Differentially Expressed Genes Table', 'Normalized Counts Table')),
                                                                                            shiny::checkboxInput(inputId = 'FilterBySignificance',
                                                                                                                 label = p(strong('Significantly Differentially Expressed')),
                                                                                                                 value = FALSE),
                                                                                            shiny::selectizeInput(inputId = 'FilterByDirection',
                                                                                                                  label = p(strong('Regulation'), br('(Choose One)')),
                                                                                                                  choices = c('Either', 'Up', 'Down'),
                                                                                                                  selected = 'Either',
                                                                                                                  multiple = FALSE,
                                                                                                                  options = list(maxItems = 1)),
                                                                                            shiny::selectizeInput(inputId = 'SearchByGene',
                                                                                                                  label = p(strong('Filter By Gene'), br('(Choose Any)')),
                                                                                                                  choices = c('All',
                                                                                                                              c(e13DRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character(),
                                                                                                                                p0DRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character(),
                                                                                                                                p10DRG_DEA_results |> dplyr::distinct(GeneID) |> unlist() |> as.character()) |> unique() ),
                                                                                                                  selected = 'All',
                                                                                                                  multiple = TRUE)
                                                                                          )
                                                                                        ))
                                                                        ),
                                    shiny::navbarMenu('Adult Hippocampus',
                                                      shiny::tabsetPanel('Comparison',
                                                                         shiny::tabPanel('Age',
                                                                                         shiny::sidebarLayout(
                                                                                           sidebarPanel = shiny::sidebarPanel(
                                                                                             shiny::selectizeInput(inputId = 'VizChoice',
                                                                                                                   label = p(strong('Visualization Type'), br('(Choose One)')),
                                                                                                                   choices = c('Volcano', 'Heatmap', 'MA Plot', 'Differentially Expressed Genes Table', 'Normalized Counts Table')),
                                                                                             shiny::checkboxInput(inputId = 'FilterBySignificance',
                                                                                                                  label = p(strong('Significantly Differentially Expressed')),
                                                                                                                  value = FALSE),
                                                                                             shiny::selectizeInput(inputId = 'FilterByDirection',
                                                                                                                   label = p(strong('Regulation'), br('(Choose One)')),
                                                                                                                   choices = c('Either', 'Up', 'Down'),
                                                                                                                   selected = 'Either',
                                                                                                                   multiple = FALSE,
                                                                                                                   options = list(maxItems = 1)),
                                                                                             shiny::selectizeInput(inputId = 'SearchByGene',
                                                                                                                   label = p(strong('Filter By Gene'), br('(Choose Any)')),
                                                                                                                   choices = c('All', ),
                                                                                                                   selected = 'All',
                                                                                                                   multiple = TRUE)
                                                                                           )
                                                                                         )),
                                                                         shiny::tabPanel('Genotype',
                                                                                         shiny::sidebarLayout(
                                                                                           sidebarPanel = shiny::sidebarPanel(
                                                                                             shiny::selectizeInput(inputId = 'VizChoice',
                                                                                                                   label = p(strong('Visualization Type'), br('(Choose One)')),
                                                                                                                   choices = c('Volcano', 'Heatmap', 'MA Plot', 'Differentially Expressed Genes Table', 'Normalized Counts Table')),
                                                                                             shiny::checkboxInput(inputId = 'FilterBySignificance',
                                                                                                                  label = p(strong('Significantly Differentially Expressed')),
                                                                                                                  value = FALSE),
                                                                                             shiny::selectizeInput(inputId = 'FilterByDirection',
                                                                                                                   label = p(strong('Regulation'), br('(Choose One)')),
                                                                                                                   choices = c('Either', 'Up', 'Down'),
                                                                                                                   selected = 'Either',
                                                                                                                   multiple = FALSE,
                                                                                                                   options = list(maxItems = 1)),
                                                                                             shiny::selectizeInput(inputId = 'SearchByGene',
                                                                                                                   label = p(strong('Filter By Gene'), br('(Choose Any)')),
                                                                                                                   choices = c('All', ),
                                                                                                                   selected = 'All',
                                                                                                                   multiple = TRUE)
                                                                                           )
                                                                                         ))
                                                                         )
                                                      )
                                    ),
                                    strong('Multiple Datasets'),
                                    shiny::tabPanel('Coming Soon')
                                    # shiny::tabPanel('6mo DRG and e13 DRG'),
                                    # shiny::tabPanel('6mo DRG and p0 DRG'),
                                    # shiny::tabPanel('6mo DRG and p10 DRG'),
                                    # shiny::tabPanel('6mo DRG and 5mo Hippocampus'),
                                    # shiny::tabPanel('6mo DRG and 15mo Hippocampus'),
                                    # shiny::tabPanel('e13 DRG and 5mo Hippocampus'),
                                    # shiny::tabPanel('e13 DRG and 15mo Hippocampus'),
                                    # shiny::tabPanel('p0 DRG and 5mo Hippocampus'),
                                    # shiny::tabPanel('p0 DRG and 15mo Hippocampus'),
                                    # shiny::tabPanel('p10 DRG and 5mo Hippocampus'),
                                    # shiny::tabPanel('p10 DRG and 15mo Hippocampus'),
                                    # shiny::tabPanel('DRG'),
                                    # shiny::tabPanel('All')))
  ),
  shiny::tabPanel(strong("IPMS"),
                  shiny::sidebarLayout(
                    sidebarPanel = shiny::sidebarPanel(
                      # shiny::selectizeInput(inputId = 'ProteomicsVizChoice',
                      #                       label = p(strong('Visualization Type'), br('(Choose One)')),
                      #                       choices = c('Volcano', 'Heatmap', 'MA Plot', 'Differentially Expressed Genes Table', 'Normalized Counts Table')),
                      shiny::checkboxInput(inputId = 'FilterByOverallEnrichment',
                                           label = p(strong('Avg TMEM FC \u+2265 2 Over Avg Control')),
                                           value = FALSE),

                      shiny::checkboxInput(inputId = 'FilterByMycEnrichment',
                                           label = p(strong('Avg TMEM-myc FC \u+2265 2 Over Avg Control-myc')),
                                           value = FALSE),

                      shiny::checkboxInput(inputId = 'FilterByV5Enrichment',
                                           label = p(strong('Avg TMEM-V5 FC \u+2265 2 Over Avg Control-V5')),
                                           value = FALSE),

                      shiny::checkboxInput(inputId = 'MWnorm',
                                           label = p(strong('Normalize Counts by Molecular Weight')),
                                           value = FALSE),

                      shiny::selectizeInput(inputId = 'FilterByHit',
                                            label = p(strong('Filter By Hit Candidate'), br('(Choose One)')),
                                            choices = c('Either', 'Up', 'Down'),
                                            selected = 'Either',
                                            multiple = FALSE,
                                            options = list(maxItems = 1)),

                      shiny::selectizeInput(inputId = 'FilterByContam',
                                            label = p(strong('Filter By Known Contaminant'), br('(Choose Any)')),
                                            choices = c('None', 'All', Crapome_results |> dplyr::filter(Num_Exp_Found > 35.8) |> dplyr::distinct(Mapped_Gene_Symbol) |> unlist() |> as.character() )),

                      shiny::selectizeInput(inputId = 'SearchByProtein',
                                            label = p(strong('Filter By Protein'), br('(Choose Any)')),
                                            choices = c('All', c(IPMS_results |> dplyr::distinct(HumanGeneID) |> unlist() |> as.character())),
                                            selected = 'All',
                                            multiple = TRUE),

                      shiny::selectizeInput(inputId = 'MapByGOTerm',
                                            label = p(strong('Map By GO Term'), br('(Choose One)')),
                                            choices = c('Default', Proteomics_GO_results |> dplyr::distinct(GO_Term) |> unlist() |> as.character()),
                                            selected = 'Default',
                                            multiple = FALSE,
                                            options = list(maxItems = 1))
                    )
                  ))
)

    # # Application title
    # titlePanel("Old Faithful Geyser Data"),
    #
    # # Sidebar with a slider input for number of bins
    # sidebarLayout(
    #     sidebarPanel(
    #         sliderInput("bins",
    #                     "Number of bins:",
    #                     min = 1,
    #                     max = 50,
    #                     value = 30)
    #     ),
    #
    #     # Show a plot of the generated distribution
    #     mainPanel(
    #        plotOutput("distPlot")
    #     )
    # )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {


  output$
    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application ----
shinyApp(ui = ui, server = server)
