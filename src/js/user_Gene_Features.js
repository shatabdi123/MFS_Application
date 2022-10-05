/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the Gene Expression menu
 */
function topmenu_change() {
    
    //Gene Structure selections    
    var struct_count = $('#mltislct option:selected').length;
    var struct_vals = $('#mltislct').val();

//GeneName - code to get the number of Genes in the box
    var genes = $('#GeneName').val();
    //console.log("genes = " + genes);
    var genes_arr = genes.replace(/(\r\n|\r|\n)/g, ',').replace(/,$/g, "").split(',');
    var gene_count = genes_arr.length;
    validate_genes(genes_arr);
   // console.log("gene_count = " + gene_count);
////////////////////////
    var analysis_select = document.getElementById('analysis');
    //GeneName - pass number of genes to populate_analysis_menu
    analysis_select.innerHTML = populate_analysis_menu(gene_count, struct_count, struct_vals);
    disable_button("bttn2", error_messages["analysis"]);
    add_error("bttn2", error_messages["analysis"]);
}

/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the analysis menu
 */
function analysis_change() {
    var analysis_val = $('#analysis option:selected').val();
    //console.log("analysis_val = " + analysis_val);
    if (analysis_val === "downsampled_hierarchial_scatter_plots") {
        var cluster = document.getElementById('cluster');
      //  console.log("cluster = " + cluster);
        if (cluster.value < 1 || !cluster.value) { console.log("setting cluster");
            cluster.setAttribute("value", 1);
            cluster.setAttribute("min", 1);
        }
        disable_button("bttn2", error_messages['analysis']);
        add_error("bttn2", error_messages["analysis"]);
    }
    else if (analysis_val == "none") {
        disable_button("bttn2", error_messages['analysis']);
        add_error("bttn2", error_messages["analysis"]);
    }
    else {
        //Enable the "Submit for analysis" buttons
        enable_button("bttn2", error_messages['analysis']);
        remove_error("bttn2", error_messages["analysis"]);
    }
}

$(document).ready(function() {
   $('#mltislct').on('change', topmenu_change);
   $('#GeneName').on('change', topmenu_change);
   $('#analysis').on('change', analysis_change);
   $('#cluster').on('change', cluster_change);
   topmenu_change();
});

/** Helper functions **/

/**
 * Populates menu options in the 'Choose analysis' select menu
 * GeneName -- updated to reflect the Gene count requirements
 */
function populate_analysis_menu(gene_count, struct_count, struct_vals) {
    const GENE_COUNT_THRESHOLD = 50;
    const MAX_GENE_STRUCTURES = 5;
    var options = "<option value='none' style='color:black'>Select one analysis</option>";
    
    if (struct_vals.includes("chromosome")) {
        options +="<optgroup label='No analyses are available for the \"Chromosome\" gene structure.'></optgroup>";
        return options;
    }
    
    options += "<option value='marginal_plots'>Marginal plot</option>";
    if (gene_count >= GENE_COUNT_THRESHOLD) {
        
        if (struct_count <= MAX_GENE_STRUCTURES) {
            options += "<optgroup label='Available for a max of " + MAX_GENE_STRUCTURES + " gene structures'>";
            options += "<option value='histogram'>Histogram</option>";
            options += "<option value='Count_and_distribution'>Count and distribution</option>";
            options += "<option value='Pair_plots'>Pair plots</option>";
            options += "<option value='Box_plots'>Box plots</option>";
            options += "<option value='Violin_plots'>Violin plots</option>";
            options += "<option value='Joint_plots'>Joint plots</option>";
            options += "<option value='Scatter_plots'>Scatter plots</option>";
            options += "<option value='Correlation_plots'>Correlation plots</option>";
            options += "</optgroup>";
            options += "<optgroup label='Available for any number of protein structures'>";
        }
        else {
            options += "<optgroup label='Reduce the number of gene stuctures to " + MAX_GENE_STRUCTURES + " or less to see additional analysis options'>";
        }
        options += "<option value='downsampled_dendogram_plots'>Downsampled Dendogram plot</option>";
        options += "<option value='downsampled_hierarchial_scatter_plots'>Downsampled Hierarchial scatterplot</option>";
        options += "<option value='downsampled_hierarchial_plots'>Downsampled Hierarchial heatmap</option>";
        options += "<option value='Codon_Heatmap_plots'>Heatmap</option>";
        options += "<option value='Codon_Gene_PCA_2D_samples_plots'>PCA 2D variables cluster</option>";
        options += "<option value='PCA_2D_label_plots'>PCA 2D observation cluster</option>";
        options += "<option value='PCA_2D_biplot_plots'>PCA 2D biplot cluster</option>";
        options += "<option value='Codon_Gene_PCA_3D_samples_plots'>PCA 3D variables cluster</option>";
        options += "<option value='PCA_3D_label_plots'>PCA 3D observation cluster</option>";
        options += "</optgroup>";
        
    }
    else {
        options += "<optgroup label='Additional analyses are available if 50 or more gene models are added'></optgroup>";
    }
    
    return options;
}

