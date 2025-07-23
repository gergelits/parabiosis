# get_parabiosis_const.r

get_pmc <- function( )
{ 
  CELLTYPES <- c( "Tconv", "CD8", "Treg", "B", "NK" )
  
  d_celltypes <- tibble(
    celltype = CELLTYPES,
    f.celltype = factor( celltype, levels = celltype ) )
  
  TISSUES <- c( "Adrenals", "BoneMarrow", "Brain", "FRT", "IEL", "Kidney",
                "Liver", "LN", "LPL", "Lung", "MLN", "Muscle",
                "Pancreas", "PP", "Skin", "Spleen", "WAT" )                     
  TISSUES.TREG <- TISSUES
  TISSUES.NK.B <- c( "Adrenals", "Brain", "FRT", "IEL", "Kidney",
                     "Liver", "LN", "LPL", "Lung", "MLN", "Muscle",
                     "Pancreas", "PP", "Skin", "Spleen", "WAT" )   
  TISSUES.NONLYMPH <- c( "Adrenals", "Brain", "FRT", "Kidney", "Liver", 
                         "Lung", "Muscle", "Pancreas", "Skin", "WAT" )   
  TISSUES.NODES <- c( "LN", "MLN" )
  TISSUES.SHORT.CIRCUIT <- c( "Spleen", "BoneMarrow" )
  TISSUES.GUT <- c( "PP", "IEL", "LPL" )
  TISSUES.NO.GUT <- TISSUES[ ! ( TISSUES %in% TISSUES.GUT ) ]
  TISSUES.NO.GUT.BM <- 
    TISSUES[ ! ( TISSUES %in% c( TISSUES.GUT, "BoneMarrow" ) ) ]
  
  tissue.all.ordered <- c( 
    "Blood", 
    "BoneMarrow", 
    "Spleen", "LN", "MLN", 
    "Adrenals", "Brain", "FRT", "Kidney", "Liver", 
    "Lung", "Muscle", "Pancreas", "Skin", "WAT", 
    "PP", "IEL", "LPL" )
  tissue.all.ordered.full <- c( 
    "Blood", 
    
    "Bone Marrow", 
    
    "Spleen", "Lymph Nodes", "Mesenteric Lymph Nodes", 
    
    "Adrenals", "Brain", "Female Reproductive Tract", "Kidney", "Liver", 
    "Lungs", "Muscle", "Pancreas", "Skin", "White Adipose Tissue", 
    
    "Peyer's Patches", "Intraepithelial Lymphocytes", 
    "Lamina Propria Lymphocytes" )    
  
  dTissueAllOrderedGroup <-
    tibble( tissue.all.ordered = tissue.all.ordered ) %>% 
      mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
      filter( tissue.all.ordered != "Blood" ) %>% 
      mutate( tissue = factor( tissue.all.ordered, levels = tissue.all.ordered ) ) %>%
      mutate( group = "Non-lymphoid") %>% 
      mutate( group = ifelse( tissue %in% c( 
        "BoneMarrow", "Spleen", "LN", "MLN", "PP"), "Lymphoid", group ) ) %>%  
      mutate( group = ifelse( tissue %in% c( "IEL", "LPL" ), "GALT", group ) ) %>% 
      mutate( group = factor( 
        group, levels = c( "Lymphoid", "Non-lymphoid", "GALT" ) ) )

  ftgl <- c( 
    "Bone marrow", "Lymphoid\ntissues", "Non-lymphoid\ntissues", 
    "Gut-assoc. lymphoid\ntissues" ) 
  
  dTissueAllOrderedGroup.f <-
    tibble( tissue.all = tissue.all.ordered ) %>% 
      mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
      filter( tissue.all != "Blood" ) %>% 
      mutate( f.tissue = factor( tissue.all, levels = tissue.all ) ) %>%
      mutate( tissue.group = "Non-lymphoid" ) %>% 
      mutate( tissue.group = ifelse( 
        f.tissue %in% "BoneMarrow", "BoneMarrow", tissue.group ) ) %>%
      mutate( tissue.group = ifelse( 
        f.tissue %in% c( "Spleen", "LN", "MLN" ), "Lymphoid", tissue.group ) ) %>%
      mutate( tissue.group = ifelse( 
        f.tissue %in% c( "PP", "IEL", "LPL" ), "GALT", tissue.group ) ) %>% 
      mutate( f.tissue.group = factor( 
        tissue.group, levels = 
          c( "BoneMarrow", "Lymphoid", "Non-lymphoid", "GALT" ) ) ) %>% 
      mutate( f.tissue.group.longer = 
                dplyr::recode( f.tissue.group,
                               BoneMarrow = ftgl[ 1 ],
                               Lymphoid = ftgl[ 2 ], 
                               `Non-lymphoid` = ftgl[ 3 ],
                               GALT = ftgl[ 4 ] ) ) %>% 
    mutate( f.tissue.group.longer = 
              factor( f.tissue.group.longer, levels = ftgl ) )
    
  
  # OLD:
  # tibble( tissue.all = tissue.all.ordered ) %>% 
  #   mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
  #   filter( tissue.all != "Blood" ) %>% 
  #   mutate( f.tissue = factor( tissue.all, levels = tissue.all ) ) %>%
  #   mutate( tissue.group = "Non-lymphoid" ) %>% 
  #   mutate( tissue.group = ifelse( 
  #     f.tissue %in% c( "BoneMarrow", "Spleen", "LN", "MLN", "PP" ), 
  #     "Lymphoid", tissue.group )) %>%
  #   mutate( tissue.group = ifelse( 
  #     f.tissue %in% c( "IEL", "LPL" ), "GALT", tissue.group ) ) %>% 
  #   mutate( f.tissue.group = factor( 
  #     tissue.group, levels = c( "Lymphoid", "Non-lymphoid", "GALT" ) ) ) %>% 
  #   mutate( f.tissue.group.longer = 
  #             dplyr::recode( f.tissue.group, 
  #                            Lymphoid = "Lymphoid\ntissues", 
  #                            `Non-lymphoid` = "Non-lymphoid\ntissues",
  #                            GALT = "Gut-assoc. lymphoid\ntissues" ) ) ->
  #   dTissueAllOrderedGroup.f

  dTissueAllOrderedGroupBlood.f <-
    tibble( tissue.all = tissue.all.ordered ) %>% 
      mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
      mutate( f.tissue = factor( tissue.all, levels = tissue.all ) ) %>%
      mutate( tissue.group = "Non-lymphoid" ) %>% 
      mutate( tissue.group = ifelse( 
        f.tissue %in% "BoneMarrow", "BoneMarrow", tissue.group ) ) %>%
      mutate( tissue.group = ifelse( 
        f.tissue %in% c( "Spleen", "LN", "MLN" ), "Lymphoid", tissue.group ) ) %>%
      mutate( tissue.group = ifelse( 
        f.tissue %in% c( "PP", "IEL", "LPL" ), "GALT", tissue.group ) ) %>% 
      mutate( tissue.group = ifelse( 
        f.tissue %in% c( "Blood" ), "Blood", tissue.group ) ) %>% 
      mutate( f.tissue.group = factor( 
        tissue.group, levels = 
          c( "BoneMarrow", "Lymphoid", "Non-lymphoid", "GALT", "Blood" ) ) ) %>% 
      mutate( f.tissue.group.longer = 
                dplyr::recode( f.tissue.group, 
                               BoneMarrow = "Bone marrow",
                               Lymphoid = "Lymphoid\ntissues", 
                               `Non-lymphoid` = "Non-lymphoid\ntissues",
                               GALT = "Gut-assoc. lymphoid\ntissues",
                               Blood = "Blood" ) )
    
    
  # OLD:
  # tibble( tissue.all = tissue.all.ordered ) %>% 
  #   mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
  #   mutate( f.tissue = factor( tissue.all, levels = tissue.all ) ) %>%
  #   mutate( tissue.group = "Non-lymphoid" ) %>% 
  #   mutate( tissue.group = ifelse( 
  #     f.tissue %in% c( "BoneMarrow", "Spleen", "LN", "MLN", "PP" ), 
  #     "Lymphoid", tissue.group )) %>%
  #   mutate( tissue.group = ifelse( 
  #     f.tissue %in% c( "IEL", "LPL" ), "GALT", tissue.group ) ) %>% 
  #   mutate( tissue.group = ifelse( 
  #     f.tissue %in% c( "Blood" ), "Blood", tissue.group ) ) %>% 
  #   mutate( f.tissue.group = factor( 
  #     tissue.group, levels = 
  #       c( "Lymphoid", "Non-lymphoid", "GALT", "Blood" ) ) ) %>% 
  #   mutate( f.tissue.group.longer = 
  #             dplyr::recode( f.tissue.group, 
  #                            Lymphoid = "Lymphoid\ntissues", 
  #                            `Non-lymphoid` = "Non-lymphoid\ntissues",
  #                            GALT = "Gut-assoc. lymphoid\ntissues",
  #                            Blood = "Blood" ) ) ->
  #   dTissueAllOrderedGroupBlood.f
  
  pmc <- list( CELLTYPES, d_celltypes,
               TISSUES, TISSUES.NONLYMPH, 
               TISSUES.NODES, TISSUES.SHORT.CIRCUIT,
               TISSUES.GUT, TISSUES.NO.GUT, TISSUES.NO.GUT.BM,
               TISSUES.NK.B, TISSUES.TREG,
               tissue.all.ordered, tissue.all.ordered.full, 
               dTissueAllOrderedGroup, dTissueAllOrderedGroup.f, 
               dTissueAllOrderedGroupBlood.f )
  names( pmc ) <- c( "CELLTYPES", "d_celltypes",
                     "TISSUES", "TISSUES.NONLYMPH", 
                     "TISSUES.NODES", "TISSUES.SHORT.CIRCUIT",
                     "TISSUES.GUT", "TISSUES.NO.GUT", "TISSUES.NO.GUT.BM",
                     "TISSUES.NK.B", "TISSUES.TREG",
                     "tissue.all.ordered", "tissue.all.ordered.full", 
                     "dTissueAllOrderedGroup", "dTissueAllOrderedGroup.f",
                     "dTissueAllOrderedGroupBlood.f" )
  return( pmc )
}


get_TISSUES <- function( celltype ) {
  if ( celltype %in% c( "Tconv", "CD8", "Treg" ) )
  {
    tissues <- get_pmc()$TISSUES
  } else if ( celltype %in% c( "NK", "B" ) )
  {
    tissues <- get_pmc()$TISSUES.NK.B
  } else 
  { 
    stop( "unknown celltype" )
  }
  return( tissues )
}
