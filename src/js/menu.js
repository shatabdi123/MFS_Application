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
var error_messages = new Array();
error_messages["one_or_no_labels"] = "One label or 'No Label' must be selected.";
error_messages["one_label"] = "Exactly one label must be selected.";
error_messages["analysis"] = "A valid analysis must be selected.";
error_messages["gene_validate"] = "One or more genes is not in the B73v5 format";
var errors = new Array();
errors["bttn1"] = "";
errors["bttn2"] = "";
errors["bttn3"] = "";
errors["analysis_container"] = "";



/**
 * Changes the status of the "Submit for analysis" and "Downsampled analysis" buttons based on the 
 * number of items selected in the Labels menu
 */
function label_change() { 
     console.log('label_select event');
    var count = $('#mltislct2 option:selected').length;
    var label = $('#mltislct2 option:selected').val();
    if (count > 1 || count == 0) {
        disable_button("bttn2", error_messages["one_or_no_labels"], true);
        disable_button("bttn3", error_messages["one_label"], true);
        disable_button("analysis_container", error_messages["one_label"], true);
        label_bttn2_enabled = false;
        label_bttn3_enabled = false;
        add_error("bttn2", error_messages["one_or_no_labels"]);
        add_error("bttn3", error_messages["one_label"]);
    }
    else if (label == "no_label") {
        if (bttn2_enabled) {
           enable_button("bttn2", error_messages["one_or_no_labels"], true);
        }
        label_bttn2_enabled = true;
        remove_error("bttn2", error_messages["one_or_no_labels"]);
        disable_button("bttn3", error_messages["one_label"], true);
        add_error("bttn3", error_messages["one_label"]);
    }
    else {
       //A non-valid option from a different page-specific menu could still be selected, 
       //if so we need to keep these items disabled
       if (bttn2_enabled) {
            enable_button("bttn2", error_messages["one_or_no_labels"], true);
       }
       label_bttn2_enabled = true;
       if (bttn3_enabled) {
            enable_button("bttn3", error_messages["one_label"], true);
            enable_button("analysis_container", error_messages["one_label"], true);
       }
       label_bttn3_enabled = true;
       remove_error("bttn2", error_messages["one_or_no_labels"]);
       remove_error("bttn3", error_messages["one_label"]);
    }
    topmenu_change();
   }

$(document).ready(function() {
     //Run the label_change() function whenever the label select menu changes
    const label_select = document.getElementById("mltislct2");
    if (label_select) {
        label_select.onchange = label_change;
        label_change();
    }
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
 * Adds a red error message about the element specified by the id
 * @param id - the button's ID
 * @param error - the error message
 */
function add_error(id, error) {
    console.log("adding error '"+error+"' to "+id);
    if (!errors[id].includes(error)) {
        errors[id] += error + "<br>";
    }
    console.log("add_error() errors " + id + ": " + errors[id]);
        $("#"+id+"_error").css("display", "block");
    $("#"+id+"_error").html(errors[id]);
}

/**
 * Adds a red error message about the element specified by the id
 * @param id - the button's ID
 * @param error - the error message
 */
function remove_error(id, error) {
    console.log("removing error '"+error+"' from "+id);
    errors[id] = errors[id].replace(error + "<br>", "");
    console.log("remove_error() errors " + id + ": " + errors[id]);
       // $("#"+id+"_error").css("display", "none");
    $("#"+id+"_error").html(errors[id]);
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

/**
 * Populate the gene text box with the selected example gene
 */
function example_GeneName(a_ele, text_ele) {
    var gene_id = a_ele.innerHTML;
    text_ele.value = gene_id;
    topmenu_change();
}

/**
 * Checks each gene one-by-one to verify they're in the valid B73v5 format
 */
function validate_genes(genes) {
  //console.log("validate_genes = " + genes + "; typeof genes = " + typeof(genes));
  for (var i=0; i<genes.length; i++) {
   // console.log("genes[i] = " + genes[i]);
    var split_test = genes[i].split("Zm00001eb");
    if (split_test[1] != null) {
        split_test[1] = split_test[1].replace(/\s/g, "");
    }
    console.log("split_test[1] = " + split_test[1] + "; typeof split_test = " + typeof(split_test) + "; split_test.length = " + split_test.length + "split_test[0] = " + split_test[0]);
    if (split_test[0] == genes[i]) {
        //Error invalid v5 gene model format
        disable_button("bttn2", error_messages["gene_validate"]);
        add_error("bttn2", error_messages["gene_validate"]);
        disable_button("bttn1", error_messages["gene_validate"]);
        add_error("bttn1", error_messages["gene_validate"]);
        return false;
    }
    if (isNaN(split_test[1]) && split_test[1] != "" && split_test[1] != null) {
     //    console.log(split_test[1] + " is not a number");
        //Error the v5 gene model does not have a number after the Zm00001eb prefix
        disable_button("bttn2", error_messages["gene_validate"]);
        add_error("bttn2", error_messages["gene_validate"]);
        disable_button("bttn1", error_messages["gene_validate"]);
        add_error("bttn1", error_messages["gene_validate"]);
        return false;
    }
  }
  enable_button("bttn2", error_messages["gene_validate"]);
  remove_error("bttn2", error_messages["gene_validate"]);
  enable_button("bttn1", error_messages["gene_validate"]);
  remove_error("bttn1", error_messages["gene_validate"]);
  return true;
}