/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the Gene Expression menu
 */
function topmenu_change() {
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
    console.log(bttn2_enabled);
    console.log(bttn3_enabled);
    console.log(label_bttn2_enabled);
    console.log(label_bttn3_enabled);
    var label = $('#mltislct2 option:selected').val();
    
     if (analysis_val == "categorical_bar_chart") {
        //Enable the "Submit for analysis" and "Downsampled analysis" buttons
        if (label_bttn2_enabled) {
            enable_button("bttn2", "");
        }
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
        disable_button("bttn3", "This button requires a label");
    }
}

$(document).ready(function() {
   $('#mltislct').on('change', topmenu_change);
   $('#analysis').on('change', analysis_change);
});

/** Helper functions **/

/**
 * Populates menu options in the 'Choose analysis' select menu
 */
function populate_analysis_menu(count, top_vals) {
    const MAX_CODONS_PLOTS = 5;
    var options = "<option value='none' style='color:black'>Select one analysis</option>";   
   
  if (compatible_analysis(top_vals)) {
       if (count <= MAX_CODONS_PLOTS) {
          options += "<optgroup label='Available for a max of " + MAX_CODONS_PLOTS + " protein localizations'>";
            options += "<option value='categorical_bar_chart'>Categorical Bar chart</option>";
          options += "</optgroup>";
          options += "<optgroup label='Available for any number of protein localizations'>";
        }
        else {
            options += "<optgroup label='Reduce the number of protein localizations to " + MAX_CODONS_PLOTS + " or less to see additional analysis options'>";
        }
        options += "<option value='Kmode_cluster_plots'>K-mode Bar plot</option>";
        options += "</optgroup>";
    }
    else {
        options += "<optgroup label='The analyses section only work for following Protein Localization features: Subcel 1-5, type and Localization.'>";
    }
    return options;
}

/**
 * Checks if any of the selected protein localizations are not compatible with
 * any of the analysis options
 *
 * @return true if no incompatible protein localizations are selected, false otherwise
 */
function compatible_analysis(localizations) {
    var loc_arr = localizations.toString().split(",");
    for (var i=0; i<loc_arr.length; i++) {
        switch (loc_arr[i]) {
            case "Nucleus":
                return false;
            case "Cytoplasm":
                return false;
            case "Extracellular":
                return false;
            case "Mitochondrion":
                return false;
            case "Cell_membrane":
                return false;
            case "Endoplasmic_reticulum":
                return false;
            case "Plastid":
                return false;
            case "Golgi_apparatus":
                return false;
            case "Lysosome_Vacuole":
                return false;
            case "Peroxisome":
                return false;
        }
    }
    return true;
}

