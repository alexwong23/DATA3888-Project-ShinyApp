library(shiny)
library(GEOquery)
library(R.utils)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(limma)
library(biomaRt)
library(ggbiplot)
library(factoextra)
library(caret)
library(glmnet)
library(DESeq2)
library(edgeR)
library(DEFormats)
library(ROCR)
library(pROC)
library(doParallel)
library(foreach)
library(DescTools)
library(ggthemes)
library(shinycssloaders)
library(shinydashboard)
library(class)
library(rsconnect)
library(BiocManager)
library(affy)
library(oligo)
# library(pd.mogene.2.0.st)
# library(mogene20sttranscriptcluster.db)
library(stringr)
# library(dashboardthemes)
library(survminer)
library(survival)
library(survMisc)
library(EnsDb.Hsapiens.v79)
options(repos = BiocManager::repositories())


survdt_m = readRDS("survdt_m.rds")

load("elastic.glm.v4.Rdata")
load("backup_model_v2.Rdata")
load("100_genes.rdata")

immune_elastic = readRDS("elastic_model_immune.RDS")
train_immune = readRDS("immune_gene_final_edit.RDS")
immune_status = readRDS("immune_status.rds")
edited_immune_elastic = readRDS("edited_elastic_model_immune.RDS")



ui <- dashboardPage(dashboardHeader(title = "DATA3888 Group 24"),
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem(text = "Risk Calculator", tabName = "RC", icon = icon("calculator")),
                        menuItem(text = "Information", tabName = "IF", icon = icon("book"))
                      )
                    ),
                    dashboardBody(

                      tabItems(

                      tabItem(tabName = "RC", h2("DATA3888 Group 24"), fluidRow(

                      box(title = "Patient Information", solidHeader = TRUE, status = "info", collapsible = TRUE,

                          radioButtons(inputId = "GenderCat", label = "Gender",
                                    choices = c("Male", "Female"),
                                    selected = "Male"
                                    ),
                          tags$p(" "),

                        selectInput(inputId = "AgeBracket", label = "Age Bracket (Years)",
                                    choices = c("25 - 35",
                                                "36 - 45",
                                                "46 - 55"),
                                    selected = "36 - 45"),

                        helpText("Please round the patient's age to the nearest whole number.
                                 For example, a 35.7 year old patient would be considered as 36 years of age. "),
                        tags$p(" "),

                        fileInput(inputId = "GeneData", label = "Patient RNA-Seq Data (.CSV)",
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv") ),

                        actionButton(inputId = "Button", label = "Calculate Risk"),

                        tags$p(" "),

                        tags$p("Please only upload a CSV file with two columns. Each row should consist of a gene and its respective expression value."),

                        tags$p("We can only accept ", tags$strong("ENSEMBL Stable IDs"), " or ",
                               tags$strong("official gene symbols approved by the HGNC. ")),

                        tags$p(tags$strong("Please do not press the button while the results are being calculated -
                                 doing so will restart the process.")),
                        
                        tags$p(" "),
                        tags$p(" "),
                        
                        tags$p(tags$strong("Data Sources")),
                        tags$p("Acute Rejection (AR) calculator is based on data from GSE120396, GSE120649, GSE131179 from Gene Expression Omnibus."),
                        tags$p("Survival curve (from 'Donor Mismatches' tab) is based on data provided by Dr. Germain Wong"),
                        tags$p("Immune Tolerance calculator is based on data from GSE120396, GSE120649, GSE131179 from Gene Expression Omnibus.")

                        
                      ),

                      tabBox(title = "Results", id = "tabset", height = 600,
                        tabPanel("AR Risk"
                                 ,plotOutput("ARPlot") %>%
                                   withSpinner(color="#0dc5c1", type = 8, size = 5),
                                 textOutput("PopMedAR")
                                 ),
                        tabPanel("Donor Mismatches"
                                 ,plotOutput("SurvPlot")
                                 ),
                        tabPanel("Immune Tolerance"
                                 ,plotOutput("ImmunePlot"),
                                 textOutput("PopMedIm")
                                 )
                      )

                    )
                    ),

                    tabItem(tabName = "IF", h2("Group 24's Risk Calculator"),fluidRow(
                      div(class = "col-sm-12 col-md-12 col-lg-12",
                      box(title = "Our Aim:", solidHeader = TRUE, status = "primary", collapsible = TRUE,width = '150px',
                          
                          
                          tags$article("Kidney organ transplantation is a lifesaving treatment for those people diagnosed with
                                       kidney disease. Kidney organ transplantation is desired over renal dialysis due to the prospect of a better a quality of life, with fewer health problems, reduced restrictions on diet and working lifestyle. "),
                          tags$article("However, kidney organ allocation has also posed itself as a major resource allocation problem. In short, there are a limited number of donor kidney organs
                                       which need to be matched effectively and accurately to the patient.  "),
                          tags$article("We aim to improve effective and accurate allocation decision-making for practitioners, as well as provide advice to potential transplant patients based on the combined interpretation of three components via a risk calculator: the genetic profile of the patient and the probability of acute rejection,
                                       eplet mismatches with donor kidneys, and the need for immunosuppression. "),
                          tags$article('Through this risk calculator, we hope that practitioners will be able to confidently assess the success of a kidney transplant for patients and make decisions regarding donor allocation. We also hope to concurrently educate patients about the factors affecting transplant success and see their
                                       predicted outcomes when compared to the general population.  ')),
          
                      box(title = "Eplets- Donor Mismatches:", solidHeader = TRUE, status = "info", collapsible = TRUE,width = '150px',
                          tags$article("Eplets are small arrangements of amino acid residues that are polymorphic - that is they occur in several forms on Human Leukocyte antigens (HLA) molecules. They are essential components of HLA epitopes that are recognized by the antibodies. HLA molecules are responsible for regulating human immune systems and therefore play an important part in graft rejection. Any cell that displays a foreign HLA type is seen as an antigen and triggers an immune response in the body which then leads to the production of antibodies,
                                       resulting in the rejecting of the tissue that bears those cells.  "),
                          tags$article('Herefore, epitope-HLA matching is important as it determines how well the match is between the donor organ and recipient patient. In essence, the greater the mismatches between the donor and the recipient
                          patient the lower the probability of survival is. By extension, the number of mismatches will better inform the practitioner the dose of
                                       immunosuppression medication needed for the patient. '),
                          tags$div('While we aim to present a general view of the epitope mismatching by stratifying the population by the characteristics of age and sex,
                                       the study of epitope compatibility is still in its infant stages and greater research needs to take place in order to drive scientific advancement',tags$em('(Tambur, 2018).')),
                          tags$article('There is a general consensus that other variables like BMI and race are important in epitope mismatching but, we have chosen to exclude these from our application. BMI has shown no significant correlation between graft survival (Papalia et.al, 2010) but, rather associated
                                       with several other socio-economic and demographic factors such as race.'),
                          tags$div('We also note that not all eplets will lead to the same type of immune response Dorr et al., 2018; Tambur, 2018).For example, two recipient-donor combinations may have identical eplet
                                       mismatches but can exhibit different immunogenicity schemes. '),
                          tags$hr(),
                          tags$cite('Tambur 2018 : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6128220/ ')),
                      box(title = " Acute Rejection & Graft Failure:", solidHeader = TRUE, status = "danger", collapsible = TRUE,width = '150px',
                          tags$article('This section looks at the top 70 differentially expressed genes used to determine whether the patient/ recipient is at risk of experiencing acute rejection. This result will be part of the discussion between patient and practitioner in their shared decision making,
                                       in addition to the other outcomes of this application. '),
                          tags$article('Graft rejection is defined as the event that the recipient’s body destroying the transplanted tissue or organ. Acute rejection is viewed as the development of antibodies that then lead to graft rejection.

Acute rejection of a transplanted organ occurs days, weeks or even months after the transplantation, where immune system of a person detects the foreign organ and seeks to attack it. Within our context, acute rejection is defined at 3 months. However, it should be noted that acute rejection is classified as the graft failing at any given time before or at three months post transplantation.
                                       It could be possible that the graft did or did not fail shortly after.  ')),
                      box(title = "Immunosuppression:", solidHeader = TRUE, status = "success", collapsible = TRUE,width = '150px',
                          tags$article('This section of the application serves to inform the practitioner and the patient, in their shared decision making, the genetic markers that indicate whether the patient will benefit from the minimizing or withdrawing from the immunosuppression medication.
                                       This will be taken into consideration along with the epitope mismatches and the genetic profile concerning acute rejection. '),
                          tags$article("Immunosuppressants are drugs that lowers the body’s immune response in order to aid the survival of the transplanted organ. Maintenance drugs like Tacrolimus and Mycophenolate Mofetil to name a couple, are medications used for the long term, ensuring prolonged graft survival. Tacrolimus in particular has
                                       been said to demonstrate better renal function (Morales et.al 2006).  "),
                          tags$article('While consuming immunosuppressant drugs are meant to aid the survival of the transplanted kidney organ, one of the major drawbacks is the increased chance
                          of infections especially following the surgery due to high dosages. Infection can further lead to rejection of the kidney and the failure of the graft. Therefore, the decision making that takes place prior to the organ transplantation is an important point of discussion between the practitioner and the patient when concerning
                                       the dosage of immunosuppressant medication as well as the duration. '),
                          tags$hr(),
                          tags$cite('Morales 2006 : https://jasn.asnjournals.org/content/17/12_suppl_3/S296?fbclid=IwAR1WCh5u-EvP4qM85aTKW3o1wPZGa6dvt_ed6Mpwr_ANMehfZk-YY09fdf0 '),
                          tags$br(),
                          tags$a(href="https://www.sciencedirect.com/topics/immunology-and-microbiology/immunosuppression","immunosuppression ,"),
                          
                          tags$a(href="https://www.kidney.org/atoz/content/immuno ","www.kidney.org")
                      ),


                          )
                    )
                    )
                    )
                    )
)

server <- function(input, output) {

  survival_plot <- eventReactive(input$Button,{

    gender_input <- input$GenderCat
    age_input <- input$AgeBracket


    user_survdt <<- survdt_m %>% dplyr::filter(
      Gender == gender_input,
      Age == age_input
    )

    user_kmcurve <<- survfit(
      Surv(C2daystodnDSA, C2dnDSA) ~ Eplet2MM, data = user_survdt)

    return(user_kmcurve)

  })

  legend_title = eventReactive(input$Button,{

    legend_title = paste(input$GenderCat, ":", input$AgeBracket, "y.o.")

    return(legend_title)
  })

  
  ar_prediction = eventReactive(input$Button,{
    
    mean_train<- mean(train_data_100)
    inFile <- input$GeneData
    
    #if (is.null(inFile))
    # return(NULL)
    
    #read.csv(inFile$datapath, header = input$header)
    file<- read.delim(inFile$datapath, header = TRUE)
    
    genes <- as.matrix(file[,1])
    gse<- as.numeric(file[,-1])
    gse<- data.matrix(gse)
    rownames(gse) <- genes
    counts <- gse
    
    #getting the processed counts matrix
    cpm_counts <- cpm(counts)
    lcpm <- log2(cpm_counts+1)
    data_counts<- normalizeQuantiles((lcpm))
    
    
    #checks if ensemble id needs to be transformed
    
    if (startsWith(genes[2], "E")){
      rownames(data_counts) = sapply(strsplit(rownames(data_counts), ".", fixed=T), function(x) x[1])
      symbols <- select(EnsDb.Hsapiens.v79, key= as.character(rownames(data_counts)), columns=c("SYMBOL"), keytype="GENEID")
      df.gse<- as.data.frame(data_counts)
      gse <- merge(df.gse, symbols, by.x = 0, by.y = 'GENEID', all.x = TRUE)
      gene_symbols<- gse[,3]
      data_counts <- as.matrix(gse[,2])
      rownames(data_counts) <- gene_symbols
    }
    
    chooseGenes <- which(!duplicated(rownames(data_counts)))
    gse<- data.matrix(data_counts[chooseGenes,])
    
    ind <- which(rownames(gse) %in% as.array(sig_genes_100))
    norm_factor <- mean_train/ mean(gse[ind,]) 
    df <- as.matrix(gse[ind,] * norm_factor)
    df2 <- df[order(rownames(df)), ]
    X = (t(df2))
    
    #us existing models
    if (length(ind)==100){
      predd_ar <- predict(elastic.glm, X, type="response")
    }
    
    
    #train another model
    if (length(ind) != 100){
      genes_present <- rownames(as.matrix(gse[ind,]))
      ind <- which(rownames(train_data_100) %in% genes_present)
      train = as.matrix(t(train_data_100[ind,]))
      y <- ifelse(rejection_status_396_179_649=="Yes", 1,0)
      data <- data.frame(y,train)
      risk = c()
      for (i in 1:5){
        
        train_cont <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "random", verboseIter = FALSE, returnData = FALSE)
        
        elastic_reg <- train( y~., data = data, method = "glmnet",
                              preProcess = c("center", "scale"),
                              tuneLength = 5,
                              trControl = train_cont
        )
        elastic.glm2 <- glmnet(train, y, alpha = elastic_reg$bestTune[1,1], lambda =elastic_reg$bestTune[1,2]	,  family = "binomial")
        
        predd_ar <- predict(elastic.glm2, X, type="response")
        risk[i] = predd_ar
        
      }
      
      predd_ar <- mean(risk)
    }
    
    return(predd_ar)

  })
  
  im_prediction = eventReactive(input$Button,{
    
    mean_train<- mean(as.matrix(train_immune))
    inFile <- input$GeneData
    
    #if (is.null(inFile))
    # return(NULL)
    
    #read.csv(inFile$datapath, header = input$header)
    file<- read.delim(inFile$datapath, header = TRUE)
    
    genes <- as.matrix(file[,1])
    gse<- as.numeric(file[,-1])
    gse<- data.matrix(gse)
    rownames(gse) <- genes
    counts <- gse
    
    #getting the processed counts matrix
    #cpm_counts <- cpm(counts)
    #lcpm <- log2(cpm_counts+1)
    #data_counts<- normalizeQuantiles((lcpm))
    
    data_counts <- counts
    
    
    #checks if ensemble id needs to be transformed
    
    if (startsWith(genes[2], "E")){
      rownames(data_counts) = sapply(strsplit(rownames(data_counts), ".", fixed=T), function(x) x[1])
      symbols <- select(EnsDb.Hsapiens.v79, key= as.character(rownames(data_counts)), columns=c("SYMBOL"), keytype="GENEID")
      df.gse<- as.data.frame(data_counts)
      gse <- merge(df.gse, symbols, by.x = 0, by.y = 'GENEID', all.x = TRUE)
      gene_symbols<- gse[,3]
      data_counts <- as.matrix(gse[,2])
      rownames(data_counts) <- gene_symbols
    }
    
    chooseGenes <- which(!duplicated(rownames(data_counts)))
    gse<- data.matrix(data_counts[chooseGenes,])
    
    ind <- which(rownames(gse) %in% as.array(rownames(train_immune)))
    norm_factor <- mean_train/ mean(gse[ind,]) 
    df <- as.matrix(gse[ind,] * norm_factor)
    df2 <- df[order(rownames(df)), ]
    X = (t(df2))
    
    #us existing models
    if (length(ind)==nrow(train_immune)){
      predd_im <- predict(immune_elastic, X, type="response")
    }
    
    
    #train another model
    if (length(ind) != nrow(train_immune)){
      genes_present <- rownames(as.matrix(gse[ind,]))
      ind <- which(rownames(train_immune) %in% genes_present)
      train = as.matrix(t(train_immune[ind,]))
      y <- ifelse(immune_status=="Immunotherapy", 0,1)
      data <- data.frame(y,train)
      risk = c()
      for (i in 1:5){
        
        train_cont <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "random", verboseIter = FALSE, returnData = FALSE)
        
        elastic_reg <- train( y~., data = data, method = "glmnet",
                              preProcess = c("center", "scale"),
                              tuneLength = 5,
                              trControl = train_cont
        )
        elastic.glm2 <- glmnet(train, y, alpha = elastic_reg$bestTune[1,1], lambda =elastic_reg$bestTune[1,2]	,  family = "binomial")
        
        predd_im <- predict(elastic.glm2, X, type="response")
        risk[i] = predd_im
        
      }
      
      predd_im <- mean(risk)
    }
    
    return(predd_im)
    
  })
  
  
  
  
  output$SurvPlot = renderPlot({

    ggsurv <- ggsurvplot(
      survival_plot(),                     # survfit object with calculated statistics.
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      palette = c("#E7B800", "#2E9FDF"),
      caption = "Data Provided by Dr. Germain Wong (University of Sydney) \n DSA = Donor Specific Antibody",
      ggtheme = theme_bw(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T,# colour risk table text annotations.
      risk.table.height = 0.25, # the height of the risk table
      risk.table.y.text = FALSE,# show bars instead of names in text annotations
      # in legend of risk table.
      surv.median.line = "hv",  # add the median survival pointer.
      xlab = "Time (Years)",
      ylab = "Estimated Probability",
      subtitle = legend_title(),
      title = "Estimated Probability for Class II de novo DSA Appearance",
      font.title = c(16, "bold", "darkblue"),
      font.x = c(14, "bold", "red"),
      font.y = c(14, "bold", "darkred"),
      font.tickslab = c(12, "bold"),
      font.subtitle = c(12, "bold"),
      legend = "top",
      legend.title = "Donor Mismatches",
      legend.labs = c(" < 30 Mismatches", " > 30 Mismatches")

    )

    ggsurv$plot =
      ggsurv$plot +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.subtitle = element_text(hjust = 0.5))

    ggsurv
  }, height = 600)
  
  
  output$ARPlot = renderPlot({
    
    df <- data.frame(matrix(nrow=2, ncol = 2))
    names(df) <- c("variable", "percentage")
    df$variable <- c("pop_risk", "pred_risk")
    df$percentage <- c(0.38, ar_prediction())
    
    df <- df %>% mutate(group=ifelse(percentage <0.28, "green",
                                     ifelse(percentage>=0.28 & percentage<0.48, "orange","red")),
                        label=paste0(round(percentage*100), "%"),
                        title=dplyr::recode(variable, `pop_risk`="Population Risk of Acute Rejection",
                                            `pred_risk`="Patient's Risk of Acute Rejection",
                        ))
    ggplot(df, aes(fill = group, ymax = percentage, ymin = 0, xmax = 2, xmin = 1)) +
      geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="#ece8bd") +
      geom_rect() + 
      coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2)) + ylim(c(0,2)) +
      geom_text(aes(x = 0, y = 0, label = label, colour=group), size=6.5) +
      geom_text(aes(x=0.8, y=1.45, label=title), size=4.2) + 
      facet_wrap(~title, ncol = 5) +
      theme_void() +
      scale_fill_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
      scale_colour_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      guides(fill=FALSE) +
      guides(colour=FALSE)+
      ggtitle("Risk of Acute Rejection")+
      theme(plot.title = element_text(size=15, face="bold"))+
      theme(plot.title = element_text(hjust = 0.5, colour = "dark blue")
      )
  })
  
  output$ImmunePlot = renderPlot({
    
    df <- data.frame(matrix(nrow=2, ncol = 2))
    names(df) <- c("variable", "percentage")
    df$variable <- c("pop_risk", "pred_risk")
    df$percentage <- c(0.85, im_prediction())
    
    df <- df %>% mutate(group=ifelse(percentage <0.20, "green",
                                     ifelse(percentage>=0.20 & percentage<0.70, "orange","red")),
                        label=paste0(round(percentage*100), "%"),
                        title=dplyr::recode(variable, `pop_risk`="Percentage of Population Relying on Immunosuppression",
                                            `pred_risk`="Patient's Predicted Reliance on Immunosuppression",
                        ))
    ggplot(df, aes(fill = group, ymax = percentage, ymin = 0, xmax = 2, xmin = 1)) +
      geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="#ece8bd") +
      geom_rect() + 
      coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2)) + ylim(c(0,2)) +
      geom_text(aes(x = 0, y = 0, label = label, colour=group), size=6.5) +
      geom_text(aes(x=0.8, y=1.45, label=title), size=4.2) + 
      facet_wrap(~title, ncol = 5) +
      theme_void() +
      scale_fill_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
      scale_colour_manual(values = c("red"="#C9146C", "orange"="#DA9112", "green"="#129188")) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      guides(fill=FALSE) +
      guides(colour=FALSE)+
      ggtitle("Reliance on Immunosuppression")+
      theme(plot.title = element_text(size=15, face="bold"))+
      theme(plot.title = element_text(hjust = 0.5, colour = "dark blue")
      )
  })


}
shinyApp(ui = ui, server = server)