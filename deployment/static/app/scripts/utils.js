'use strict';

jQuery('img.svg').each(function(){
    var $img = jQuery(this);
    var imgID = $img.attr('id');
    var imgClass = $img.attr('class');
    var imgURL = $img.attr('src');

    jQuery.get(imgURL, function(data) {
        // Get the SVG tag, ignore the rest
        var $svg = jQuery(data).find('svg');

        // Add replaced image's ID to the new SVG
        if(typeof imgID !== 'undefined') {
            $svg = $svg.attr('id', imgID);
        }
        // Add replaced image's classes to the new SVG
        if(typeof imgClass !== 'undefined') {
            $svg = $svg.attr('class', imgClass+' replaced-svg');
        }

        // Remove any invalid XML tags as per http://validator.w3.org
        $svg = $svg.removeAttr('xmlns:a');

        // Replace image with new SVG
        $img.replaceWith($svg);

    }, 'xml');

});

/* Small function to apply a tick to steps which have been completed */
function applyTicks(current_step) {
    if(current_step == "intro" || current_step == "add") {
        $('span.glyphicon.add').addClass('hidden');
        $('span.glyphicon.map').addClass('hidden');
        $('span.glyphicon.validate').addClass('hidden');
    }
    else if(current_step == "map") {
        $('span.glyphicon.add').removeClass('hidden');
        $('span.glyphicon.map').addClass('hidden');
        $('span.glyphicon.validate').addClass('hidden');
    }
    else if (current_step == "validate") {
        $('span.glyphicon.add').removeClass('hidden');
        $('span.glyphicon.map').removeClass('hidden');
        $('span.glyphicon.validate').addClass('hidden');
    }
    else if(current_step == "finish") {
        $('span.glyphicon.add').removeClass('hidden');
        $('span.glyphicon.map').removeClass('hidden');
        $('span.glyphicon.validate').removeClass('hidden');
    }
}

function splitInput(in_str) {
    var data = [];
    var splitted = in_str.split('\n');
    if (splitted.length == 1){
        //try delimiting on commas
        splitted = in_str.split(' ');
    }

    //return some json
    jQuery.each(splitted, function(index, value) {

        //work out if it's SMILES or InChi (for now)
        var input_type = smilesOrInchi(value);
        data.push({ type: input_type || "unknown" , input_str: value });
    });

    return data;
}

function smilesOrInchi(in_str) {
    //if it starts with InChi, it's inchi
    //otherwise SMILES (for now)
    var type_str = ""
    //is it inchi?
    if (in_str.trim().match(/^((InChI=)?[^J][0-9BCOHNSOPrIFla+\-\(\)\\\/,pqbtmsih]{6,})$/ig)) {
        type_str = "InChi";
    }
    //is it inchikey?
    else if (27===in_str.length && '-'===in_str[14] && '-'===in_str[25] && !!in_str.match(/^([0-9A-Z\-]+)$/)) {
        type_str = "InChi Key"
    }
    //is it a chembl id?
    else if (in_str.trim().match(/^((CHEMBL)?[0-9])$/)) {
        type_str = "ChEMBL ID"
    }
    else {
        type_str = "SMILES"
    }
    return type_str;
}