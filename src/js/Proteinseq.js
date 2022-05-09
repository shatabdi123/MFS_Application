console.log('Proteinseq.js loaded');

/**
 * Enables or disables the "Submit for analysis" and "Downsampled analysis" buttons based on the value of the DNA Composition and Autocorrelation menu
 */
function topmenu_change() {
    var comp_value = $('#selUser option:selected').val();
    console.log("comp_value = " + comp_value);
    var analysis_select = document.getElementById('selUser1');
    analysis_select.innerHTML = populate_analysis_menu(comp_value);
    $("#selUser1").select2("destroy");
    $("#selUser1").select2();
    if (comp_value == "DC" || comp_value == "TC") {
        //Enable the "Submit for analysis" and "Downsampled analysis" buttons
        if (label_bttn2_enabled) {
            enable_button("bttn2", "");
        }
        if (label_bttn3_enabled) {
            enable_button("bttn3", "");
        }
        enable_button("selUser1", "");
        enable_button("analysis_container", "");
    }
    else {
        //Enable the "Submit for analysis" and "Downsampled analysis" buttons
        disable_button("bttn2", "This button can only be used if 'dinucleic' or 'trinucleic' are selected in the 'DNA Composition and Autocorrelation' menu");
        disable_button("bttn3", "This button can only be used if 'dinucleic' or 'trinucleic' are selected in the 'DNA Composition and Autocorrelation' menu");
        disable_button("selUser1", "This button can only be used if 'dinucleic' or 'trinucleic' are selected in the 'DNA Composition and Autocorrelation' menu");
        disable_button("analysis_container", "Can only be used if 'dinucleic' or 'trinucleic' are selected in the 'DNA Composition and Autocorrelation' menu");
    }
}
$(document).ready(function() {
    //Run the comp_change() function whenever the 'DNA Composition and Autocorrelation' menu changes
    var comp_select = document.getElementById("selUser");
    console.log("comp_select = " + comp_select);
   // comp_select.onchange = comp_change;
   $('#selUser').on('change', topmenu_change);
   
   
});

/** Helper functions **/

/**
 * Populates menu options in the 'Choose analysis' select menu
 */
function populate_analysis_menu(type) {
    var options = "<optgroup label='Graphical analysis'>";
    if (type == "DC") {
        options += "<option value='count_dipeptide_box_plots'>Frequency Dipeptide box plot</option>";
    }
    else if (type == "TC") {
        options += "<option value='count_tripeptide_box_plots'>Frequency Tripeptide box plot</option>";
    }
    else {
        options += "<option value='none_"+type+"' style='color:black' selected>Analyses are for dipeptide or tripeptide compositions</option>";
    }
    options += "</optgroup>";
    return options;
}