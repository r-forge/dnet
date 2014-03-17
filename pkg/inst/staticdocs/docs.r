library(staticdocs)
library(grid)
list(
    readme = "",
    
    index = list(
        sd_section("Analysis functions",
            "These functions are used for analysis",
            c(
                'eCal',
                'eView',
                'eWrite',
                'dPvalAggregate',
                'dNetInduce',
                'dBUMfit',
                'dBUMscore',
                'dNetFind',
                'dNetPipeline',
                'dNetConfidence',
                'dNetReorder',
                'dRWR',
                'dContrast',
                'dCommSignif',
                'dSVDSignif',
                'dFDRscore'
            )
        ),
        sd_section("Visualisation functions",
            "These functions are used for visualisation",
            c(
                'visRunES',
                'visNet',
                'visNetMul',
                'visNetReorder',
                'visNetArc',
                'visNetCircle',
                'visHeatmap'
            )
        ),
        sd_section("Built-in databases",
            "",
            c(
                "org.Mm.eg",
                "org.Mm.egDO", 
                "org.Mm.egGOBP",
                "org.Mm.egGOCC",
                "org.Mm.egGOMF",
                "org.Mm.egHPMI", 
                "org.Mm.egHPON",
                "org.Mm.egHPPA",
                "org.Mm.egMP",
                "org.Mm.egPS",
                "org.Hs.string",
                "org.Mm.string",
                "org.At.string"
            )
        ),
        sd_section("Built-in datasets",
            "",
            c(
                "Hiratani_TableS1",
                "CLL"
            )
        )
    ),
    
    if(0){
    icons = list(  
        eCal = sd_icon({
          textGrob("Common", rot = 45, gp = gpar(cex = 1))
        }),
        visRunES = sd_icon({
          textGrob("Hot", rot = 45, gp = gpar(cex = 1.2))
        }),
        eView = sd_icon(inherit = "eCal"),
        visNet = sd_icon(inherit = "visRunES")
    )
    }
)
