
/**
 * Hack function to deal with two top menus on the TFbindingSite page.
 * This function is only called when the label changes in menu.js
 *
 * Calls the analysis_change() function to check the new value of label
 */
function topmenu_change() {
    analysis_change();
}

/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the Gene Expression menu
 */
function topmenu_singleslct_change() {
console.log("singleslct change");
    $("#mltislct option:selected").prop("selected", false);
    $("#mltislct").multiselect("destroy");
    $('#mltislct').multiselect({
        maxHeight: 200,
        buttonWidth:'310px',
        enableFiltering: true,
        enableCaseInsensitiveFiltering: true,
        includeSelectAllOption: false,
        filterPlaceholder:'Search Here..'
    });
    var analysis_select = document.getElementById('analysis');
    analysis_select.innerHTML = populate_analysis_menu(0, "");
    disable_button("bttn2", "An analysis must be selected");
    disable_button("bttn3", "An analysis must be selected");
}


/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the Gene Expression menu
 */
function topmenu_multislct_change() {
    var top_singleslct = document.getElementById('top_singleslct');
    top_singleslct.value = "none";
    
    var top_count = $('#mltislct option:selected').length;
    var top_vals = $('#mltislct').val();
    console.log("top_count= " + top_count);
    console.log("top_vals= " + top_vals);
    var analysis_select = document.getElementById('analysis');
    analysis_select.innerHTML = populate_analysis_menu(top_count, top_vals);
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
   $('#top_singleslct').on('change', topmenu_singleslct_change);
   $('#mltislct').on('change', topmenu_multislct_change);
   $('#analysis').on('change', analysis_change);
   $('#cluster').on('change', cluster_change);
});

/** Helper functions **/

/**
 * Populates menu options in the 'Choose analysis' select menu
 */
function populate_analysis_menu(count, top_vals) {
    const MAX_TOP_COUNT = 5;
    var options = "<option value='none' style='color:black'>Select one analysis</option>";
   
    if (count <= MAX_TOP_COUNT) {
      options += "<optgroup label='Available for a max of " + MAX_TOP_COUNT + " gene structures'>";
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
        options += "<optgroup label='Reduce the number of gene stuctures to " + MAX_TOP_COUNT + " or less to see additional analysis options'>";
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

