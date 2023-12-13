# function to create dataframe with overlaps (small)
create_stacked_dataframe <- function(region_1, region_2, region_3, region_4) {
  stacked_dataframe <- data.frame(
    "brain_region" = c(
      rep(region_1, 7), rep(region_2, 7),
      rep(region_3, 7), rep(region_4, 7)
    ),
    "Analysis" = rep(c(
      "DGE only", "DTE only", "DTU only", "DGE and DTE",
      "DGE and DTU", "DTE and DTU", "all"
    ), 4),
    "n_genes" = c(
      # dge only
      length(setdiff(
        cerebellum_dge_genes,
        union(
          intersect(cerebellum_dge_genes, cerebellum_dte_genes),
          intersect(cerebellum_dge_genes, cerebellum_dtu_genes)
        )
      )),
      # dte only
      length(setdiff(
        cerebellum_dte_genes,
        union(
          intersect(cerebellum_dge_genes, cerebellum_dte_genes),
          intersect(cerebellum_dte_genes, cerebellum_dtu_genes)
        )
      )),
      # dtu only
      length(setdiff(
        cerebellum_dtu_genes,
        union(
          intersect(cerebellum_dge_genes, cerebellum_dtu_genes),
          intersect(cerebellum_dte_genes, cerebellum_dtu_genes)
        )
      )),
      # dge and dte
      length(intersect(
        cerebellum_dge_genes, cerebellum_dte_genes
      )) -
        length(intersect(
          intersect(cerebellum_dge_genes, cerebellum_dte_genes),
          cerebellum_dtu_genes
        )),
      # dge and dtu
      length(intersect(cerebellum_dge_genes, cerebellum_dtu_genes)) -
        length(intersect(
          intersect(cerebellum_dge_genes, cerebellum_dte_genes),
          cerebellum_dtu_genes
        )),
      # dte and dtu
      length(intersect(cerebellum_dtu_genes, cerebellum_dte_genes)) -
        length(intersect(
          intersect(cerebellum_dge_genes, cerebellum_dte_genes),
          cerebellum_dtu_genes
        )),
      # all
      length(intersect(
        intersect(cerebellum_dge_genes, cerebellum_dte_genes),
        cerebellum_dtu_genes
      )),
      
      # dge only
      length(setdiff(cortex_dge_genes, union(intersect(
        cortex_dge_genes,
        cortex_dte_genes
      ), intersect(cortex_dge_genes, cortex_dtu_genes)))),
      # dte only
      length(setdiff(cortex_dte_genes, union(intersect(
        cortex_dge_genes,
        cortex_dte_genes
      ), intersect(cortex_dte_genes, cortex_dtu_genes)))),
      # dtu only
      length(setdiff(cortex_dtu_genes, union(intersect(
        cortex_dge_genes,
        cortex_dtu_genes
      ), intersect(cortex_dte_genes, cortex_dtu_genes)))),
      # dge and dte
      length(intersect(cortex_dge_genes, cortex_dte_genes)) - length(intersect(
        intersect(cortex_dge_genes, cortex_dte_genes),
        cortex_dtu_genes
      )),
      # dge and dtu
      length(intersect(cortex_dge_genes, cortex_dtu_genes)) - length(intersect(
        intersect(cortex_dge_genes, cortex_dte_genes),
        cortex_dtu_genes
      )),
      # dte and dtu
      length(intersect(cortex_dtu_genes, cortex_dte_genes)) - length(intersect(
        intersect(cortex_dge_genes, cortex_dte_genes),
        cortex_dtu_genes
      )),
      # all
      length(intersect(
        intersect(cortex_dge_genes, cortex_dte_genes),
        cortex_dtu_genes
      )),
      
      # dge only
      length(setdiff(hippocampus_dge_genes, union(
        intersect(hippocampus_dge_genes, hippocampus_dte_genes),
        intersect(hippocampus_dge_genes, hippocampus_dtu_genes)
      ))),
      # dte only
      length(setdiff(hippocampus_dte_genes, union(
        intersect(hippocampus_dge_genes, hippocampus_dte_genes),
        intersect(hippocampus_dte_genes, hippocampus_dtu_genes)
      ))),
      # dtu only
      length(setdiff(hippocampus_dtu_genes, union(
        intersect(hippocampus_dge_genes, hippocampus_dtu_genes),
        intersect(hippocampus_dte_genes, hippocampus_dtu_genes)
      ))),
      # dge and dte
      length(intersect(hippocampus_dge_genes, hippocampus_dte_genes)) -
        length(intersect(
          intersect(hippocampus_dge_genes, hippocampus_dte_genes),
          hippocampus_dtu_genes
        )),
      # dge and dtu
      length(intersect(hippocampus_dge_genes, hippocampus_dtu_genes)) -
        length(intersect(
          intersect(hippocampus_dge_genes, hippocampus_dte_genes),
          hippocampus_dtu_genes
        )),
      # dte and dtu
      length(intersect(hippocampus_dtu_genes, hippocampus_dte_genes)) -
        length(intersect(
          intersect(hippocampus_dge_genes, hippocampus_dte_genes),
          hippocampus_dtu_genes
        )),
      # all
      length(intersect(
        intersect(hippocampus_dge_genes, hippocampus_dte_genes),
        hippocampus_dtu_genes
      )),
      
      # dge only
      length(setdiff(striatum_dge_genes, union(
        intersect(
          striatum_dge_genes, striatum_dte_genes
        ),
        intersect(
          striatum_dge_genes, striatum_dtu_genes
        )
      ))),
      # dte only
      length(setdiff(striatum_dte_genes, union(
        intersect(
          striatum_dge_genes, striatum_dte_genes
        ),
        intersect(
          striatum_dte_genes, striatum_dtu_genes
        )
      ))),
      # dtu only
      length(setdiff(striatum_dtu_genes, union(
        intersect(
          striatum_dge_genes, striatum_dtu_genes
        ),
        intersect(
          striatum_dte_genes, striatum_dtu_genes
        )
      ))),
      # dge and dte
      length(intersect(striatum_dge_genes, striatum_dte_genes)) -
        length(intersect(
          intersect(striatum_dge_genes, striatum_dte_genes),
          striatum_dtu_genes
        )),
      # dge and dtu
      length(intersect(striatum_dge_genes, striatum_dtu_genes)) -
        length(intersect(
          intersect(striatum_dge_genes, striatum_dte_genes),
          striatum_dtu_genes
        )),
      # dte and dtu
      length(intersect(striatum_dtu_genes, striatum_dte_genes)) -
        length(intersect(
          intersect(striatum_dge_genes, striatum_dte_genes),
          striatum_dtu_genes
        )),
      # all
      length(intersect(
        intersect(striatum_dge_genes, striatum_dte_genes),
        striatum_dtu_genes
      ))
    )
  )
  
  return(stacked_dataframe)
}

# function to create a long dataframe for stacked barplots
create_stacked_dataframe_paired <- function(comparison_1, comparison_2,
                                            comparison_3, comparison_4,
                                            comparison_5, comparison_6){
  stacked_barplot_paired <- data.frame(
    "brain_regions" = c(
      rep(comparison_1, 7),
      rep(comparison_2, 7),
      rep(comparison_3, 7),
      rep(comparison_4, 7),
      rep(comparison_5, 7),
      rep(comparison_6, 7)
    ),
    "Analysis" = rep(c(
      "DGE only", "DTE only", "DTU only", "DGE and DTE",
      "DGE and DTU", "DTE and DTU", "all"
    ), 6),
    "n_genes" = c(
      # dge only
      length(setdiff(
        cerebellum_cortex_dge_genes,
        union(
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dte_genes
          ),
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dtu_genes
          )
        )
      )),
      # dte only
      length(setdiff(
        cerebellum_cortex_dte_genes,
        union(
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dte_genes
          ),
          intersect(
            cerebellum_cortex_dte_genes,
            cerebellum_cortex_dtu_genes
          )
        )
      )),
      # dtu only
      length(setdiff(
        cerebellum_cortex_dtu_genes,
        union(
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dtu_genes
          ),
          intersect(
            cerebellum_cortex_dte_genes,
            cerebellum_cortex_dtu_genes
          )
        )
      )),
      # dge and dte
      length(intersect(
        cerebellum_cortex_dge_genes,
        cerebellum_cortex_dte_genes
      )) -
        length(intersect(
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dte_genes
          ),
          cerebellum_cortex_dtu_genes
        )),
      # dge and dtu
      length(intersect(
        cerebellum_cortex_dge_genes,
        cerebellum_cortex_dtu_genes
      )) -
        length(intersect(
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dte_genes
          ),
          cerebellum_cortex_dtu_genes
        )),
      # dte and dtu
      length(intersect(
        cerebellum_cortex_dtu_genes,
        cerebellum_cortex_dte_genes
      )) -
        length(intersect(
          intersect(
            cerebellum_cortex_dge_genes,
            cerebellum_cortex_dte_genes
          ),
          cerebellum_cortex_dtu_genes
        )),
      # all
      length(intersect(
        intersect(
          cerebellum_cortex_dge_genes,
          cerebellum_cortex_dte_genes
        ),
        cerebellum_cortex_dtu_genes
      )),
      
      # dge only
      length(setdiff(
        cerebellum_hippocampus_dge_genes,
        union(
          intersect(
            cerebellum_hippocampus_dge_genes,
            cerebellum_hippocampus_dte_genes
          ),
          intersect(
            cerebellum_hippocampus_dge_genes,
            cerebellum_hippocampus_dtu_genes
          )
        )
      )),
      # dte only
      length(setdiff(
        cerebellum_hippocampus_dte_genes,
        union(
          intersect(
            cerebellum_hippocampus_dge_genes,
            cerebellum_hippocampus_dte_genes
          ),
          intersect(
            cerebellum_hippocampus_dte_genes,
            cerebellum_hippocampus_dtu_genes
          )
        )
      )),
      # dtu only
      length(setdiff(
        cerebellum_hippocampus_dtu_genes,
        union(
          intersect(
            cerebellum_hippocampus_dge_genes,
            cerebellum_hippocampus_dtu_genes
          ),
          intersect(
            cerebellum_hippocampus_dte_genes,
            cerebellum_hippocampus_dtu_genes
          )
        )
      )),
      # dge and dte
      length(intersect(
        cerebellum_hippocampus_dge_genes,
        cerebellum_hippocampus_dte_genes
      )) - length(intersect(
        intersect(
          cerebellum_hippocampus_dge_genes,
          cerebellum_hippocampus_dte_genes
        ),
        cerebellum_hippocampus_dtu_genes
      )),
      # dge and dtu
      length(intersect(
        cerebellum_hippocampus_dge_genes,
        cerebellum_hippocampus_dtu_genes
      )) - length(intersect(
        intersect(
          cerebellum_hippocampus_dge_genes,
          cerebellum_hippocampus_dte_genes
        ),
        cerebellum_hippocampus_dtu_genes
      )),
      # dte and dtu
      length(intersect(
        cerebellum_hippocampus_dtu_genes,
        cerebellum_hippocampus_dte_genes
      )) - length(intersect(
        intersect(
          cerebellum_hippocampus_dge_genes,
          cerebellum_hippocampus_dte_genes
        ),
        cerebellum_hippocampus_dtu_genes
      )),
      # all
      length(intersect(
        intersect(
          cerebellum_hippocampus_dge_genes,
          cerebellum_hippocampus_dte_genes
        ),
        cerebellum_hippocampus_dtu_genes
      )),
      
      # dge only
      length(setdiff(
        cerebellum_striatum_dge_genes,
        union(
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dte_genes
          ),
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dtu_genes
          )
        )
      )),
      # dte only
      length(setdiff(
        cerebellum_striatum_dte_genes,
        union(
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dte_genes
          ),
          intersect(
            cerebellum_striatum_dte_genes,
            cerebellum_striatum_dtu_genes
          )
        )
      )),
      # dtu only
      length(setdiff(
        cerebellum_striatum_dtu_genes,
        union(
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dtu_genes
          ),
          intersect(
            cerebellum_striatum_dte_genes,
            cerebellum_striatum_dtu_genes
          )
        )
      )),
      # dge and dte
      length(intersect(
        cerebellum_striatum_dge_genes,
        cerebellum_striatum_dte_genes
      )) -
        length(intersect(
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dte_genes
          ),
          cerebellum_striatum_dtu_genes
        )),
      # dge and dtu
      length(intersect(
        cerebellum_striatum_dge_genes,
        cerebellum_striatum_dtu_genes
      )) -
        length(intersect(
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dte_genes
          ),
          cerebellum_striatum_dtu_genes
        )),
      # dte and dtu
      length(intersect(
        cerebellum_striatum_dtu_genes,
        cerebellum_striatum_dte_genes
      )) -
        length(intersect(
          intersect(
            cerebellum_striatum_dge_genes,
            cerebellum_striatum_dte_genes
          ),
          cerebellum_striatum_dtu_genes
        )),
      # all
      length(intersect(
        intersect(
          cerebellum_striatum_dge_genes,
          cerebellum_striatum_dte_genes
        ),
        cerebellum_striatum_dtu_genes
      )),
      
      # dge only
      length(setdiff(
        cortex_hippocampus_dge_genes,
        union(
          intersect(
            cortex_hippocampus_dge_genes,
            cortex_hippocampus_dte_genes
          ),
          intersect(
            cortex_hippocampus_dge_genes,
            cortex_hippocampus_dtu_genes
          )
        )
      )),
      # dte only
      length(setdiff(
        cortex_hippocampus_dte_genes,
        union(
          intersect(
            cortex_hippocampus_dge_genes,
            cortex_hippocampus_dte_genes
          ),
          intersect(
            cortex_hippocampus_dte_genes,
            cortex_hippocampus_dtu_genes
          )
        )
      )),
      # dtu only
      length(setdiff(
        cortex_hippocampus_dtu_genes,
        union(
          intersect(
            cortex_hippocampus_dge_genes,
            cortex_hippocampus_dtu_genes
          ),
          intersect(
            cortex_hippocampus_dte_genes,
            cortex_hippocampus_dtu_genes
          )
        )
      )),
      # dge and dte
      length(intersect(
        cortex_hippocampus_dge_genes,
        cortex_hippocampus_dte_genes
      )) -
        length(intersect(intersect(
          cortex_hippocampus_dge_genes,
          cortex_hippocampus_dte_genes
        ), cortex_hippocampus_dtu_genes)),
      # dge and dtu
      length(intersect(
        cortex_hippocampus_dge_genes,
        cortex_hippocampus_dtu_genes
      )) - length(intersect(intersect(
        cortex_hippocampus_dge_genes,
        cortex_hippocampus_dte_genes
      ), cortex_hippocampus_dtu_genes)),
      # dte and dtu
      length(intersect(
        cortex_hippocampus_dtu_genes,
        cortex_hippocampus_dte_genes
      )) - length(intersect(intersect(
        cortex_hippocampus_dge_genes,
        cortex_hippocampus_dte_genes
      ), cortex_hippocampus_dtu_genes)),
      # all
      length(intersect(
        intersect(
          cortex_hippocampus_dge_genes,
          cortex_hippocampus_dte_genes
        ),
        cortex_hippocampus_dtu_genes
      )),
      
      # dge only
      length(setdiff(
        cortex_striatum_dge_genes,
        union(
          intersect(
            cortex_striatum_dge_genes,
            cortex_striatum_dte_genes
          ),
          intersect(cortex_striatum_dge_genes, cortex_striatum_dtu_genes)
        )
      )),
      # dte only
      length(setdiff(
        cortex_striatum_dte_genes,
        union(
          intersect(
            cortex_striatum_dge_genes,
            cortex_striatum_dte_genes
          ),
          intersect(
            cortex_striatum_dte_genes,
            cortex_striatum_dtu_genes
          )
        )
      )),
      # dtu only
      length(setdiff(
        cortex_striatum_dtu_genes,
        union(
          intersect(
            cortex_striatum_dge_genes,
            cortex_striatum_dtu_genes
          ),
          intersect(
            cortex_striatum_dte_genes,
            cortex_striatum_dtu_genes
          )
        )
      )),
      # dge and dte
      length(intersect(cortex_striatum_dge_genes, cortex_striatum_dte_genes)) -
        length(intersect(
          intersect(
            cortex_striatum_dge_genes,
            cortex_striatum_dte_genes
          ),
          cortex_striatum_dtu_genes
        )),
      # dge and dtu
      length(intersect(cortex_striatum_dge_genes, cortex_striatum_dtu_genes)) -
        length(intersect(
          intersect(
            cortex_striatum_dge_genes,
            cortex_striatum_dte_genes
          ),
          cortex_striatum_dtu_genes
        )),
      # dte and dtu
      length(intersect(cortex_striatum_dtu_genes, cortex_striatum_dte_genes)) -
        length(intersect(
          intersect(
            cortex_striatum_dge_genes,
            cortex_striatum_dte_genes
          ),
          cortex_striatum_dtu_genes
        )),
      # all
      length(intersect(
        intersect(
          cortex_striatum_dge_genes,
          cortex_striatum_dte_genes
        ),
        cortex_striatum_dtu_genes
      )),
      
      # dge only
      length(setdiff(
        striatum_hippocampus_dge_genes,
        union(
          intersect(
            striatum_hippocampus_dge_genes,
            striatum_hippocampus_dte_genes
          ),
          intersect(
            striatum_hippocampus_dge_genes,
            hippocampus_striatum_dtu_genes
          )
        )
      )),
      # dte only
      length(setdiff(
        hippocampus_striatum_dte_genes,
        union(
          intersect(
            striatum_hippocampus_dge_genes,
            striatum_hippocampus_dte_genes
          ),
          intersect(
            hippocampus_striatum_dte_genes,
            hippocampus_striatum_dtu_genes
          )
        )
      )),
      # dtu only
      length(setdiff(
        hippocampus_striatum_dtu_genes,
        union(
          intersect(
            striatum_hippocampus_dge_genes,
            hippocampus_striatum_dtu_genes
          ),
          intersect(
            striatum_hippocampus_dte_genes,
            hippocampus_striatum_dtu_genes
          )
        )
      )),
      # dge and dte
      length(intersect(
        striatum_hippocampus_dge_genes,
        striatum_hippocampus_dte_genes
      )) -
        length(intersect(
          intersect(
            striatum_hippocampus_dge_genes,
            striatum_hippocampus_dte_genes
          ),
          hippocampus_striatum_dtu_genes
        )),
      # dge and dtu
      length(intersect(
        striatum_hippocampus_dge_genes,
        hippocampus_striatum_dtu_genes
      )) -
        length(intersect(
          intersect(
            striatum_hippocampus_dge_genes,
            striatum_hippocampus_dte_genes
          ),
          hippocampus_striatum_dtu_genes
        )),
      # dte and dtu
      length(intersect(
        hippocampus_striatum_dtu_genes,
        striatum_hippocampus_dte_genes
      )) -
        length(intersect(
          intersect(
            striatum_hippocampus_dge_genes,
            striatum_hippocampus_dte_genes
          ),
          hippocampus_striatum_dtu_genes
        )),
      # all
      length(intersect(
        intersect(
          striatum_hippocampus_dge_genes,
          striatum_hippocampus_dte_genes
        ),
        hippocampus_striatum_dtu_genes
      ))
    )
  )
  
  return(stacked_barplot_paired)
}