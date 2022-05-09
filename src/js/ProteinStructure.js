/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the Gene Expression menu
 */
function topmenu_change() {
    var exp_count = $('#mltislct option:selected').length;
    var exp_vals = $('#mltislct').val();
    console.log("exp_count= " + exp_count);
    console.log("exp_vals= " + exp_vals);
    var analysis_select = document.getElementById('analysis');
    analysis_select.innerHTML = populate_analysis_menu(exp_count, exp_vals);
    disable_button("bttn2", "An analysis must be selected");
    disable_button("bttn3", "An analysis must be selected");
}

/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the analysis menu
 */
function analysis_change() {
    var analysis_val = $('#analysis option:selected').val();
    console.log("analysis_val = " + analysis_val);
    var label = $('#mltislct2 option:selected').val();
     if (isLimitedPlotAnalysis(analysis_val)){
        //Enable the "Submit for analysis" and "Downsampled analysis" buttons
        if (label_bttn2_enabled) {
            enable_button("bttn2", "");
        }
        if (label_bttn3_enabled) {
            enable_button("bttn3", "");
        }

    }
    else if (analysis_val === "down_hierarchial_scatter_plots") {
        var cluster = document.getElementById('cluster');
        console.log("cluster = " + cluster);
        if (cluster.value < 1 || !cluster.value) { console.log("setting cluster");
            cluster.setAttribute("value", 1);
            cluster.setAttribute("min", 1);
        }
        disable_button("bttn2", "This button can't be used with the currently selected analysis");
        if (label_bttn3_enabled) {
            enable_button("bttn3", "");
        }
    }
    else if (analysis_val == "none") {
        disable_button("bttn2", "An analysis must be selected");
        disable_button("bttn3", "An analysis must be selected");
    }
    else {
        //Enable the "Submit for analysis" and "Downsampled analysis" buttons
        disable_button("bttn2", "This button can't be used with the currently selected analysis");
        if (label_bttn3_enabled) {
            enable_button("bttn3", "");
        }
    }
    
    if (label == "no_label") {
        //enable_button("bttn2", "");
        disable_button("bttn3", "This button requires a label");
    }
}

$(document).ready(function() {
   $('#mltislct').on('change', topmenu_change);
   $('#analysis').on('change', analysis_change);
   $('#cluster').on('change', cluster_change);
});

/** Helper functions **/

/**
 * Populates menu options in the 'Choose analysis' select menu
 */
function populate_analysis_menu(count, exp_vals) {
    const MAX_CODONS_PLOTS = 5;
    var options = "<option value='none' style='color:black'>Select one analysis</option>";
    console.log('populate_analysis menu, exp_vals = ' + exp_vals);
     if (exp_vals == "SignalP") {
        //If Signal Peptides is the only selection     
        options += "<optgroup label='Available for \"Signal Peptides\" protein structure only'>"
        options += "<option value='categorical_bar_chart'>Categorical Bar chart</option></optgroup>";
        return options;
     }
     else if (exp_vals.includes("SignalP")) {
     console.log('SignalP is included');
        //If Signal Peptides is one of multiple selections
        options += "<optgroup label='Signal Peptides can not be selected with other protein structures'>"
        options += "</optgroup>";
        return options;
       
     }

    if (count <= MAX_CODONS_PLOTS) {
      options += "<optgroup label='Available for a max of " + MAX_CODONS_PLOTS + " protein structures'>";
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
        options += "<optgroup label='Reduce the number of protein stuctures to " + MAX_CODONS_PLOTS + " or less to see additional analysis options'>";
    }
    options += "<option value='down_dendogram_plots'>Downsampled Dendogram plot</option>";
    options += "<option value='down_hierarchial_scatter_plots'>Downsampled Hierarchial scatterplot</option>";
    options += "<option value='down_hierarchial_plots'>Downsampled Hierarchial heatmap</option>";
    
    options += "<option value='Count_Heatmap_plots'>Heatmap</option>";
    options += "<option value='Count_Gene_PCA_2D_samples_plots'>PCA 2D variables cluster</option>";
    options += "<option value='PCA_2D_label_plots'>PCA 2D observation cluster</option>";
    options += "<option value='PCA_2D_biplot_plots'>PCA 2D biplot cluster</option>";
    options += "<option value='Count_Gene_PCA_3D_samples_plots'>PCA 3D variables cluster</option>";
    options += "<option value='PCA_3D_label_plots'>PCA 3D observation cluster</option>";
    options += "</optgroup>";
    return options;
}

