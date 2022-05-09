/**
 * menu.js
 * Purpose: Holds generic functions used by multiple MFS pages
 */
console.log('menu.js loaded');

/**
 * Global variables to check the status of the buttons from the page-specific javascript functions
 */
var bttn2_enabled = true;
var bttn3_enabled = true;
var label_bttn2_enabled = true;
var label_bttn3_enabled = true;



/**
 * Changes the status of the "Submit for analysis" and "Downsampled analysis" buttons based on the 
 * number of items selected in the Labels menu
 */
function label_change() { 
     console.log('label_select event');
    var count = $('#mltislct2 option:selected').length;
    var label = $('#mltislct2 option:selected').val();
    console.log("num labels = " + count);
    if (count > 1 || count == 0) {
        disable_button("bttn2", "This feature can only be used when one label is selected.", true);
        disable_button("bttn3", "This feature can only be used when one label is selected.", true);
        disable_button("analysis_container", "This feature can only be used when one label is selected.", true);
        label_bttn2_enabled = false;
        label_bttn3_enabled = false;
    }
    else if (label == "no_label") {
        if (bttn2_enabled) {
           enable_button("bttn2", "", true);
        }
        label_bttn2_enabled = true;
        disable_button("bttn3", "This feature cannot be used with no labels.", true);
    }
    else {
       //A non-valid option from a different page-specific menu could still be selected, 
       //if so we need to keep these items disabled
       if (bttn2_enabled) {
            enable_button("bttn2", "", true);
       }
       label_bttn2_enabled = true;
       if (bttn3_enabled) {
            enable_button("bttn3", "", true);
            enable_button("analysis_container", "", true);
       }
       label_bttn3_enabled = true;
    }
    topmenu_change();
   }

$(document).ready(function() {
     //Run the label_change() function whenever the label select menu changes
    const label_select = document.getElementById("mltislct2");
    label_select.onchange = label_change;
    label_change();
});
    

/***** Helper Functions ********/

/**
 * Enables a button with an optional mouse-over message
 * @param id - the button's ID
 * @param msg - an optional mousover message
 * @param label_menu - indicates whether this function was called from the label menu changing or a different menu
 */
 
function enable_button(id, msg="", label_menu=false) {
    $("#"+id).prop("disabled",false);
    $("#"+id).css("cursor", "pointer");
    $("#"+id).prop("title", msg);
    $("#"+id).css("filter", "brightness(100%)");
    
    if (!label_menu) {
        if (id == "bttn2") {
            bttn2_enabled = true;
        }
        else if (id == "bttn3") {
            bttn3_enabled = true;
        }
    }
}

/**
 * Disables a button with an optional mouse-over message
 * @param id - the button's ID
 * @param msg - an optional mousover message
 * @param label_menu - indicates whether this function was called from the label menu changing or a different menu
 */
function disable_button(id, msg="", label_menu=false) {
    $("#"+id).prop("disabled",true);
    $("#"+id).css("cursor","not-allowed");
    $("#"+id).css("filter", "brightness(70%)");
    $("#"+id).prop("title", msg);
    
     if (!label_menu) {
        if (id == "bttn2") {
            bttn2_enabled = false;
        }
        else if (id == "bttn3") {
            bttn3_enabled = false;
        }
     }
}

/**
 * Checks if the provided analysis is one of the ones that are limited to 5 Codon selections.
 * Written for the Codon page, but also used on the GeneExp
 *
 * @param analysis - the analysis type to be checked
 * @return true if the analysis is limited by 5 Codon selections, false otherwise
 */
function isLimitedPlotAnalysis(analysis) {
  var limited_analyses = new Array("histogram", "Count_and_distribution", "Pair_plots", "Box_plots", "Violin_plots", "Joint_plots", "Scatter_plots", "Correlation_plots");
  for (var i=0; i<limited_analyses.length; i++) {
    if (analysis === limited_analyses[i]) {
        return true;
    }
  }
  return false;
}

/**
 * Forces a min and max value for the cluster input
 */
function cluster_change() {
    var cluster = document.getElementById('cluster');
    const MAX_CLUSTER = 10;
    const MIN_CLUSTER = 1;
    if (cluster.value > MAX_CLUSTER) {
        cluster.value = MAX_CLUSTER;
    }
    else if (cluster.value < MIN_CLUSTER) {
        cluster.value = MIN_CLUSTER;
    }
}