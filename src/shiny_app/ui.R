# shiny app and paper overview tab
welcome_about <- tabPanel(
  title = "Welcome and About",
  br(),
  br(),
  br(),
  br(),
  br(),
  div(
    h3("Welcome!"),
    p(
      class = "welcome-about-p",
      "Welcome to the Shiny application created by Emma Jones in the",
      a("Lasseigne Lab", href = "https://www.lasseigne.org/"),
      "for visualizing long-read mouse brain RNA-sequencing data!"
    ),
    p(
      class = "welcome-about-p",
      "To address differences in splicing across brain regions and sexes, we used long-read ONT RNA sequencing to sequence 40 mouse brain cDNA libraries from 10 mice and calculated differential transcript usage. We found that there is strong evidence of differential transcript usage across brain regions as well as differential expression at the gene and transcript level. We found that the brain region with the most differential expression and transcript usage is the cerebellum, potentially driven by differences in cell type composition. Overall, our findings suggest there is much differential splicing across brain regions and to a lesser extent, within brain regions across sexes. These differences in splicing could explain sex differences in prevalence and prognosis of various neurological and psychiatric disorders."
    ),
    tags$figure(
      align = "center",
      tags$img(
      src = "graphical_abstract.png",
      width = "825px",
      height = "225px"
    ),
    tags$figcaption(
      "Graphical overview of our data generation and analysis",
      "If you are interested in checking out our work, feel free to check out our accompanying manuscript",
      a("here", href = "placeholder.com"),
      ".")
    ),
    p(
      class = "welcome-about-p",
      "This web application allows users to query their own genes of interest in our generated dataset and produce plots to explore their own hypotheses. Happy learning!"
    ),
    hr(),
    h3("What is differential transcript usage?"),
    p(
      class = "welcome-about-p",
      "Differential transcript usage (DTU) measures changes in the proportions of transcripts expressed that a gene may have across groups."
    ),
    tags$figure(
      align = "center",
      tags$img(
      src = "DTU_example.png",
      width = "600px",
      height = "300px"
    ),
    tags$figcaption("An example of classical differential gene expression on the left, and differential transcript usage on the right.")
    ),
    p(
      class = "welcome-about-p",
      "In the figure above, there is an example on the left of differential gene expression across conditions, but no differential transcript usage since the transcript isoforms are still expressed in the same proportions. On the right, you can see that the transcript isoforms are expressed in different proportions, showing differential transcript usage without differential gene expression."
    ),
    p(
      class = "welcome-about-p",
      "The isoform switch plots displayed are generated using the",
      a("IsoformSwitchAnalyzeR",
        href = "https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html"),
      "package created and maintained by Kristoffer Vitting-Seerup."
    ),
    hr(),
    h3("Contact Information"),
    p(
      class = "welcome-about-p",
      "For any issues or inquiries about this website or its related publication, please contact Emma Jones (efjones@uab.edu) or Brittany Lasseigne (bnp0001@uab.edu)."
    ),
    hr(),
    div(
      class = "branding",
      p("All figures created with BioRender.com"),
      p(
        "Copyright 2023 by the",
        tags$a(href = "https://www.lasseigne.org/", "Lasseigne Lab")
      )
    )
  )
)

# gene expression heatmap tab
gene_exp_heatmap <- tabPanel(
  title = "Custom Gene Expression Heatmap",
  br(),
  br(),
  br(),
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      p(
        "Interested in gene-level expression (counts per million) of your genes of interest?"
      ),
      p(
        "Use this tool to make a heatmap of you favorite genes (max recommended genes is around 30), and see how they cluster colored by tissue sample."
      ),
      selectizeInput(
        inputId = "gene_ids",
        label = "Input Gene Symbols or ENSEMBL IDs Here",
        choices = NULL,
        multiple = TRUE
      ),
      downloadButton("download_heatmap", "Download Raw Data"),
      downloadButton("download_heatmap_image", "Download Image", icon = icon("camera")),
      br(),
      br(),
      p(
        "Note: be mindful of expression values, heatmap color may be skewed when a very lowly expressed gene is plotted with more highly expressed genes."
      )
    ),
    mainPanel(imageOutput("heatmap_image"))
  )
)

# single brain region vs all others switchplot tab
single_region_plot <- tabPanel(
  title = "Compare Single Brain Region to Others",
  br(),
  br(),
  br(),
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      p(
        "Are you interested in comparing differential transcript usage across a single brain region compared to other regions?"
      ),
      p(
        "If you are, select a gene name and a brain region to focus on. Alternatively, you can download the processed data to plot yourself!"
      ),
      selectizeInput(
        inputId = "gene_name1",
        label = "Input Gene Symbol Here",
        choices = NULL
      ),
      selectInput(
        inputId = "brain_region_1",
        label = "Select Brain Region",
        choices = list(
          "Cerebellum" = "cerebellum",
          "Cortex" = "cortex",
          "Hippocampus" = "hippocampus",
          "Striatum" = "striatum"
        ),
        selected = "cerebellum"
      ),
      downloadButton("download_data_1", "Download Raw Data"),
      downloadButton("download_image_1", "Download Image", icon = icon("camera")),
      br(),
      br(),
      p(
        "Note: if you do not see your gene of interest, it is likely it did not have at least 2 transcripts measured in our data or was not expressed in both samples and/or brain regions."
      )
    ),
    mainPanel(imageOutput("switchplot_1_image"))
  )
)

# pairwise brain region switchplot tab
pairwise_region_plot <- tabPanel(
  title = "Compare Two Brain Regions",
  br(),
  br(),
  br(),
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      p(
        "Are you interested in comparing differential transcript usage pairwise across two specific brain regions?"
      ),
      p(
        "If you are, select a gene name and brain regions to focus on. Alternatively, you can download the processed data to plot yourself!"
      ),
      selectizeInput(
        inputId = "gene_name2",
        label = "Input Gene Symbol Here",
        choices = NULL
      ),
      selectInput(
        inputId = "condition1",
        label = "Select Brain Region 1",
        choices = list(
          "Cerebellum" = "cerebellum",
          "Cortex" = "cortex",
          "Hippocampus" = "hippocampus",
          "Striatum" = "striatum"
        ),
        selected = "cerebellum"
      ),
      selectInput(
        inputId = "condition2",
        label = "Select Brain Region 2",
        choices = list(
          "Cerebellum" = "cerebellum",
          "Cortex" = "cortex",
          "Hippocampus" = "hippocampus",
          "Striatum" = "striatum"
        ),
        selected = "cortex"
      ),
      downloadButton("download_data_2", "Download Raw Data"),
      downloadButton("download_image_2", "Download Image", icon = icon("camera")),
      br(),
      br(),
      p(
        "Note: if you do not see your gene of interest, it is likely it did not have at least 2 transcripts measured in our data or was not expressed in both samples and/or brain regions."
      )
    ),
    mainPanel(imageOutput("switchplot_2_image"))
  )
)

# within brain region across sexes switchplot tab
sex_region_plot <- tabPanel(
  title = "Compare Sex Within a Brain Region",
  br(),
  br(),
  br(),
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      p(
        "Are you interested in comparing differential transcript usage within a single brain region across sexes?"
      ),
      p(
        "If you are, select a gene name and a brain region to focus on. Alternatively, you can download the processed data to plot yourself!"
      ),
      selectizeInput(
        inputId = "gene_name3",
        label = "Input Gene Symbol Here",
        choices = NULL
      ),
      selectInput(
        inputId = "brain_region_3",
        label = "Select Brain Region",
        choices = list(
          "Cerebellum" = "cerebellum",
          "Cortex" = "cortex",
          "Hippocampus" = "hippocampus",
          "Striatum" = "striatum"
        ),
        selected = "cerebellum"
      ),
      downloadButton("download_data_3", "Download Raw Data"),
      downloadButton("download_image_3", "Download Image", icon = icon("camera")),
      br(),
      br(),
      p(
        "Note: if you do not see your gene of interest, it is likely it did not have at least 2 transcripts measured in our data or was not expressed in both samples and/or brain regions."
      )
    ),
    
    mainPanel(imageOutput("switchplot_3_image"))
  )
)

ui <- navbarPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    tags$script(type = "text/javascript", src = "code.js")
  ),
  title = "Visualizing Isoform Switches in Wild Type Mouse Brain",
  position = "fixed-top",
  welcome_about,
  gene_exp_heatmap,
  single_region_plot,
  pairwise_region_plot,
  sex_region_plot
)