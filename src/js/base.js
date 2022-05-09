

/**
 * Submit search form
 */
function search_submit(search_text) {
   document.location = 'https://www.google.com/search?hl=en&ie=ISO-8859-1&btnG=Google+Search&q=site%3Ahttps://mfs.maizegdb.org+' + search_text;
}


/**
* Check if enter key was pressed, and submit search form if so
*/
function check_key(e) {
    if (e.key === "Enter") {
        search_submit(document.getElementById('search_box').value);
    }
}



