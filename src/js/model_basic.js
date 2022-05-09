

/**
 * Forces a min and max value for the cluster input
 */
function validate_range(ele, num_decimals) {
    console.log("ele = " + ele);
    var val = Number.parseFloat(ele.value);
    var min = Number.parseFloat(ele.min);
    var max = Number.parseFloat(ele.max);
 
    if (val > max) {
        ele.value = max.toFixed(num_decimals);
    }
    else if (val < min) {
        ele.value = min.toFixed(num_decimals);
    }
    else if (Number.isInteger(val)) {
        console.log('integer found, parsing as float');
        ele.value = val.toFixed(num_decimals);
    }
}

function validate_sequence(ele, seq_type) {
    console.log("ele = " + ele);
    var arr = ele.value.replace(/\s/g, "").toUpperCase().split("");
    var submit_btn = document.getElementById("seq_submit");
    if (seq_type == 'dna') {
       var error_div = document.getElementById('dna_seq_error');
       for (var i=0; i<arr.length; i++) {
            var c = arr[i];
            if (c != "A" && c != "G" && c != "C" && c != "T" && c != "U" && c != "N") {
                //Dna sequence is invalid
                submit_btn.disabled = true;
                error_div.style.display = "inline";
                document.getElementById('dna_char_error').innerHTML = c;
                return false;
            }
       }
       submit_btn.disabled = false;
       error_div.style.display = "none";
       return true;
    }
    else if (seq_type == 'prot') {
       var error_div = document.getElementById('prot_seq_error');
       for (var i=0; i<arr.length; i++) {
            var c = arr[i];
            if (c != "A" && c != "C" && c != "D" && c != "E" && c != "F" && c != "G" &&
                c != "H" && c != "I" && c != "J" && c != "K" && c != "L" && c != "M" &&
                c != "N" && c != "O" && c != "P" && c != "Q" && c != "R" && c != "S" &&
                c != "T" && c != "V" && c != "W" && c != "Y" && c != "X" && c != "*" &&
                c != "." && c != "Z" && c != "B") {
                
                //Prot sequence is invalid
                submit_btn.disabled = true; console.log('invalid prot seq');
                error_div.style.display = "inline";
                document.getElementById('prot_char_error').innerHTML = c;
                return false;
            }
       } console.log('valid prot seq');
       submit_btn.disabled = false;
       error_div.style.display = "none";
       return true;
    }

}

function example_dnaseq(ele) {
    
var seq = `CGATCGATCGAGATATCGTGTCAACGTGCCGGCCGGGGCGTGGGAAGATGATGAACCTATCGGCTGCCGC
CAACGGCCGCGACGAGTTCCCCCCCTACGTCGTGCCGTCCAACGCGGCCGCTCCGCCCCCTTCCCTGCTC
CCAACCATGGAGCAGCAGCAGGAGAGCAGCATCCACAGGGAGCATCATCAGCTGCTGGGCTACAACCTCG
AGGCCAACTCGCTGGCCCTCCTGCCCCCGTCCAACGCCGCCGCCGCCCACCACCACACCACCTTCGCCGG
CGGCCACAGCCCCCACGACATCCTCCACTTCTACACACCTCCTCCTTCCGCCGCCTCGCACTACCTCGCC
GCCGCCGCCGGCAACCCCTACAGCCACTTAGTCTCCGCGCCCGGGACCACCTTCCACCAGACCTCGTCGT
CCTACTACCCGCCCGCGGCGGCGGCGCAGGCCGCGCCCGAGTACTACTTCCCCACCCTCGTCAGCTCCGC
CGAGGAGAACATGGCCAGCTTCGCCGCCACGCAGCTCGGCCTCAACCTCGGCTACCGCACCTACTTCCCG
CCCAGAGGAGGCTACACGTACGGCCACCACCCGCCGCGCTGCCAGGCCGAGGGCTGCAAGGCCGACCTCT
CCAGCGCCAAGCGATACCACCGTCGCCACAAGGTGTGCGAGCACCACTCCAAGGCGCCCGTCGTCGTCAC
CGCCGGTGGACTGCATCAGAGGTTCTGCCAGCAGTGCAGCAGATTCCATCTGCTGGATGAGTTCGACGAT
GCTAAGAAGAGCTGCAGGAAACGCCTAGCGGACCACAACCGCCGCCGCCGGAAGTCAAAGCCATCGGATG
CTGATGCCGGAGACAAGAAAAGGGCACATGCGAACAAAGCAGCTGCAGCTAAAGACAAAGCAGAGAGTAG
CAGCAAGAACATGGATATCGGAGATGGGTTAGGCGCACAGATACTGGGAAGTGCACTCTTGTCCAAGGAA
CAAGATCAAACCATGGATCTTGGAGAGGTGGTGAAAGAAGCAGTGGATCCCAAGGGGAAGGCATCAATGC
AACAGCATTACGGCTTCCCCTTCCATTCGTCGTCAGCAGGATCTTGCTTCCCCCAGACCCAAGCCGTCTC
CAGTGATACCACATCCAATATAGGTCAAGTGCAAGAGCCAAGCTTAGGGTTCCACCATCAGCACCACCAA
CACAGCAACATCTTGCAGCTCGGCCAGGCTATGTTTGATCTCGACTTCGATCACTAGTCAATATGTGATG
CACATGCACCCTCTCTTTCTCACCCCACCCCACCCCTCCCCTCCCCTCTATCTTGTTTGTGCGCGTAATC
CGAATATTTTTCCTTTTTAAATTATCTGTGTTAATTACTGTAACGTGGACATAATAATGATAGTCTATGC
TTGCCATGCAA`;
    ele.value = seq;
}

function example_protseq(ele) {
    var seq = `MFYSHQLLARKAPLGQIWMAATLHSKINRKRLDKLDIIKICEEILNPSVPMALRLSGILMGGVVIVYERK
VKLLYTDVSRLLTEINEAWRIKPVTDPTVLPKGKTQAKYEAVTLPEINMVVEQPMFFSEPDGAKFRRMGL
EDLDEQYVQVNLDDDDFSHADDRHQAKAVNITLVDNFESGLAETDLFNHFERFDIADDETTVNITPDEYP
QVPSTLIPSPPRQEDIPQQEEPYYAAPSPVHGEPQQGGPEDQEEQKMKQPPKASKRKARWEVPRVIMDNN
QMMIPGNIYQTWLKDASSLVSKRRKLNSNFNFIRSTKISDLMHIPPVALISHDNLFSELCYPKPLMQLWK
DCTEVKSTKASSGGQRSSSQEPQPKNSPPQAGGEYEMETGGLPMDLTDGIEKLRANMSAKYDRAYNILHS
DHSVTPGSPGLSRRSASSSGGSGSAFIQLDPEVQLPSGSGRSKRGQHSSARSLGNLDTVEEDFPLEQEVR
DFKMRRLSDYVPTPDLLEETEPTQTPYERRSNPMDKITETIQSHLKLHFDTPGVPQSESLSHLAHGMTKA
RAARLFYQIAVLATCDYIKVTQLERKGDELYGDILISRGLKM`;
    ele.value = seq;
}


/*
$(document).ready(function() {

});
*/

/** Helper functions **/


