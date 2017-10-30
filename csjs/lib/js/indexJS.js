/*
 ## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
 ## Centre National de la Recherche Scientifique>

 ## 1. This program is free software: you can redistribute it and/or modify
 ## it under the terms of the GNU Affero General Public License as published
 ## by the Free Software Foundation version 3 of the  License and under the
 ## terms of article 2 below.
 ## 2. This program is distributed in the hope that it will be useful, but
 ## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 ## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General
 ## Public License for more details.
 ## You should have received a copy of the GNU Affero General Public License
 ## along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ## 3. Communication to the public by any means, in particular in the form of
 ## a scientific paper, a poster, a slideshow, an internet page, or a patent,
 ## of a result obtained directly or indirectly by running this program must
 ## cite the following paper :
 ##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
 ##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
 ##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
 ##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
 ##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
 ## -------------------------------------------------------------------------

 ## Authors (alphabetically): Jacob L., Jaillard M., Lima L.
 */

//these globals are used in function longColumnRenderer()
var maxLengthColumnRenderer = 20;
var pathToLib = "components/lib/"


//create the preview of a component
function buildComponentPreview(componentId, preview, annotationsHOT) {
    //add the preview and hide it
    $("#showItemsDiv").append(preview);
    $("#table_comp_"+componentId.toString()).hide();

    //get the annotation div
    var annotationDiv = $('#annot_comp_'+componentId.toString());

    //if there is no annotation, then it is easy
    if (annotationsHOT.length==0) {
        annotationDiv.html("<b>No annotations found.</b>")
    }
    else {
      //now, we gotta build the HOT containing the annotation
      //populate the table for the graph annotation
      var annotationTableSettings = {
          data: annotationsHOT,
          columns: [
              {renderer: longColumnRenderer},
              {type: 'text'},
              {type: 'text'}
          ],
          colHeaders: [
              'Annotation',
              '# nodes',
              'E-value'
          ],
          copyColsLimit: 1000000,
          copyRowsLimitNumber: 1000000,
          readOnly: true,
          wordWrap: false,
          stretchH: 'all',
          allowInsertColumn: false,
          allowInsertRow: false,
          allowRemoveColumn: false,
          allowRemoveRow: false,
          autoColumnSize: {useHeaders: true},
          autoWrapCol: true,
          autoWrapRow: true,
          manualColumnResize: true,
          columnSorting: true,
          sortIndicator: true,
          height: 200
      };
      
      
      var annotationTableContainer = document.getElementById('annot_comp_'+componentId.toString());
      annotationTable = new Handsontable(annotationTableContainer, annotationTableSettings);
      annotationTable.sort(1, false);
    }
}

//create all the components preview
function buildAllComponents() {
    var idPreviewAnnOfAllTables = alasql('SELECT id, preview, annHOT FROM Components');

    idPreviewAnnOfAllTables.forEach(function(idPreviewAnn){
        console.log(idPreviewAnn)
        //buildComponentPreview
    });
}


function buildObjectsHTML(objects) {
    var html="";
    var nb=0;
    objects.forEach(function(object) {
        html += object.preview + "\n";
        nb += 1;
    });
    return "<h3>Done! Found " + nb + " components.</h3><br/><br/>" + html;
}

function showObjects() {
    var objectDiv = $('#showItemsDiv');
    objectDiv.html('Working...');

    //get the sql
    var sqlAfterSelect = ""
    if ($('#builder').queryBuilder('getSQL', false).sql.length == 0) {
        sqlAfterSelect = ' FROM Components ORDER BY ' + $('#sortTerm').val() + ' ' + $('#asc_desc').val();
    }else {
        sqlAfterSelect = ' FROM Components WHERE ' + $('#builder').queryBuilder('getSQL', false).sql + ' ORDER BY ' + $('#sortTerm').val() + ' ' + $('#asc_desc').val();
    }

    //retrieve the selected objects
    var objects = alasql('SELECT preview' + sqlAfterSelect);
    objectDiv.html(buildObjectsHTML(objects));

    //execute the js to fill in the annotations
    var annotationsJS = alasql('SELECT annScript' + sqlAfterSelect);

    annotationsJS.forEach(function(annotationJS){
        eval(annotationJS.annScript)
    });
}


//load the data in the data variable into the AlaSQL DB
function createComponentsDB(data) {
  alasql('CREATE TABLE Components');
  alasql.tables.Components.data = data;
}



//javascript function to collapse and uncolappse the stat figures
$(function(){
    $('#arguments').on('hide.bs.collapse', function () {
        $('#argumentsButton').html('<span class="glyphicon glyphicon-collapse-down"></span> Show arguments used to produce these results');
    })
    $('#arguments').on('show.bs.collapse', function () {
        $('#argumentsButton').html('<span class="glyphicon glyphicon-collapse-up"></span> Hide arguments used to produce these results');
    })
})

$(function(){
    $('#filters').on('hide.bs.collapse', function () {
        $('#filtersButton').html('<span class="glyphicon glyphicon-collapse-down"></span> Show filters');
    })
    $('#filters').on('show.bs.collapse', function () {
        $('#filtersButton').html('<span class="glyphicon glyphicon-collapse-up"></span> Hide filters');
    })
})

$(function(){
    $('#stats').on('hide.bs.collapse', function () {
        $('#statsButton').html('<span class="glyphicon glyphicon-collapse-down"></span> Show figures on lineage effect');
    })
    $('#stats').on('show.bs.collapse', function () {
        $('#statsButton').html('<span class="glyphicon glyphicon-collapse-up"></span> Hide figures on lineage effect');
    })
})
