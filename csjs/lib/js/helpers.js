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

//****************************************************************
//declare some global variables
var cy; //cytoscape object
var nodeTableData = []; //node table data
var table; //the table
var widthTable = parseInt(0.8*$(window).width(), 10);
var heightTable= parseInt(0.2*$(window).height(), 10);
var timeout = null

//some variables to deal with user interaction
var dialogGetFASTANodes;
var dialogInstructions;
var cytoscapeExportDialog;
var cytoscapeDesktopGraph = null;
var colors = ['red', 'blue', 'green', 'yellow', 'fuchsia', 'brown', 'lime', 'aqua', 'Aquamarine', 'BlueViolet', 'CadetBlue', 'DarkBlue', 'DarkCyan', 'DarkGoldenRod', 'DarkGray', 'DarkGreen', 'DarkMagenta', 'DarkSlateBlue', 'DarkSeaGreen', 'DarkSalmon', 'DarkViolet', 'DarkTurquoise'];
//declare some global variables
//****************************************************************



//*******************************************************
//FUNCTIONS CONCERNING CLIPBOARD
function copyTextToClipboard(text) {
    var textArea = document.createElement("textarea");

    //
    // *** This styling is an extra step which is likely not required. ***
    //
    // Why is it here? To ensure:
    // 1. the element is able to have focus and selection.
    // 2. if element was to flash render it has minimal visual impact.
    // 3. less flakyness with selection and copying which **might** occur if
    //    the textarea element is not visible.
    //
    // The likelihood is the element won't even render, not even a flash,
    // so some of these are just precautions. However in IE the element
    // is visible whilst the popup box asking the user for permission for
    // the web page to copy to the clipboard.
    //

    // Place in top-left corner of screen regardless of scroll position.
    textArea.style.position = 'fixed';
    textArea.style.top = 0;
    textArea.style.left = 0;

    // Ensure it has a small width and height. Setting to 1px / 1em
    // doesn't work as this gives a negative w/h on some browsers.
    textArea.style.width = '2em';
    textArea.style.height = '2em';

    // We don't need padding, reducing the size if it does flash render.
    textArea.style.padding = 0;

    // Clean up any borders.
    textArea.style.border = 'none';
    textArea.style.outline = 'none';
    textArea.style.boxShadow = 'none';

    // Avoid flash of white box if rendered for any reason.
    textArea.style.background = 'transparent';


    textArea.value = text;

    document.body.appendChild(textArea);

    textArea.select();

    try {
        var successful = document.execCommand('copy');
        var msg = successful ? 'successful' : 'unsuccessful';
    } catch (err) {
        alert('Failed to copy to clipboard...');
    }

    document.body.removeChild(textArea);
}

function getFastaOfSelectedNodes() {
    var selectedNodes = cy.$('node:selected');
    var str="";
    for (i = 0; i < selectedNodes.length; i++)
        str += ">" + selectedNodes[i].id() + "\n" + selectedNodes[i].data('name') + "\n";
    return str;
}

function copyFastaToClipboard() {
    var str = getFastaOfSelectedNodes();
    copyTextToClipboard(str);
}

//FUNCTIONS CONCERNING CLIPBOARD
//*******************************************************



//*******************************************************
//FUNCTIONS CONCERNING COLORS
function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}
//FUNCTIONS CONCERNING COLORS
//*******************************************************


//************************************************************
//FUNCTIONS OF THE CONTEXT MENU
function selectAllNodes() {
    cy.nodes().select();
};

function unselectAllNodes() {
    cy.nodes().unselect();
};
//FUNCTIONS OF THE CONTEXT MENU
//************************************************************

//************************************************************
//FUNCTIONS OF DRAWING
function drawGradient(){
    var c = document.getElementById("gradientCanvas");
    var ctx = c.getContext("2d");

    var grd = ctx.createLinearGradient(0, 0, 150, 0);
    grd.addColorStop(0, "blue");
    grd.addColorStop(1, "red");

    ctx.fillStyle = grd;
    ctx.fillRect(0, 0, 150, 50);
}


function drawAlleles() {
    var c = document.getElementById("alleleCanvas");
    var ctx = c.getContext("2d");

    ctx.beginPath();
    ctx.arc(50, 50, 50, 0, 2 * Math.PI, false);
    ctx.fillStyle = 'white';
    ctx.fill();

    ctx.arc(175, 50, 5, 0, 2 * Math.PI, false);
    ctx.fillStyle = 'white';
    ctx.fill();
}
//FUNCTIONS OF DRAWING
//************************************************************

//************************************************************
//FUNCTIONS OF SELECT/UNSELECT

var selectedNodeEqualsToNodeTableDataNodeFunction = function selectedNodeEqualsToNodeTableDataNode(selectedNode, nodeTableDataNode) {
    return nodeTableDataNode[0] == selectedNode.id()
}

var nodeTableDataNodeEqualsToSelectedNodeFunction = function nodeTableDataNodeEqualsToSelectedNode(nodeTableDataNode, selectedNode) {
    return nodeTableDataNode[0] == selectedNode.id()
}

function alreadyPresentInCollection(element, collection, testFunction) {
    for (var i = 0; i < collection.length; i++) {
        if(testFunction(element, collection[i]))
            return true;
    }
    return false;
}

function fillTable() {
    if (timeout!=null)
        clearTimeout( timeout );
    timeout = setTimeout(function(){
        var selectedNodes = cy.$('node:selected');

        //flag which nodes are already present in nodeTableData
        var alreadyPresentInNodeTableData= []
        selectedNodes.forEach(function(node){
            alreadyPresentInNodeTableData.push(alreadyPresentInCollection(node, nodeTableData, selectedNodeEqualsToNodeTableDataNodeFunction))
        })

        //add the nodes to the table if they are not yet there
        for (var i = 0; i < selectedNodes.length; i++) {
            node = selectedNodes[i]
            if (alreadyPresentInNodeTableData[i] == false) {
                nodeTableData.push([node.id(),
                    node.data('total'),
                    node.data('pheno0'),
                    node.data('pheno1'),
                    //TODO: NA temporarily removed
                    //node.data('NA'),
                    node.data('tags').toString(),
                    node.data('significant'),
                    node.data('qValue'),
                    node.data('weight'),
                    node.data('waldStatistic'),
                    node.data('sequenceLength'),
                    node.data('name')
                ])
            }
        }

        //remove from nodeTableData the nodes that are not anymore selected
        nodeTableDataTemp = []
        nodeTableData.forEach(function(node){
            if (alreadyPresentInCollection(node, selectedNodes, nodeTableDataNodeEqualsToSelectedNodeFunction) == false) {
                //do nothing
            }else {
                nodeTableDataTemp.push(node);
            }
        })

        //clear nodeTableData
        nodeTableData.length = 0;

        //fill nodeTableData with correct nodes
        nodeTableDataTemp.forEach(function(node){
            nodeTableData.push(node);
        })

        //render the table
        table.render();
    }, 100); // may have to adjust this val
}



function selectNodesFromATag (option, cy, DBGWAS_graph_tag2nodes) {
  var tag = $(option).val();

  if (tag=="clear") {
    unselectAllNodes();
    cy.fit()
  }else {
    var nodesToHighlight=[];

    cy.nodes().forEach(function (node) {
        if (DBGWAS_graph_tag2nodes[tag].includes(node.id()))
            nodesToHighlight.push(node)
    })
    
    unselectAllNodes();
    cy.collection(nodesToHighlight).select();
    cy.fit(cy.collection(nodesToHighlight))
  }
}
//FUNCTIONS OF SELECT/UNSELECT
//************************************************************



//************************************************************
//FUNCTIONS OF EXPORTING
//function to export the graph to Cytoscape Desktop
function makeFile (text, file, fileType) {
    var data = new Blob([text], {type: fileType});

    // If we are replacing a previously generated file we need to
    // manually revoke the object URL to avoid memory leaks.
    if (file !== null) {
        window.URL.revokeObjectURL(file);
    }

    file = window.URL.createObjectURL(data);

    return file;
}
//FUNCTIONS OF EXPORTING
//************************************************************




//************************************************************
//FUNCTIONS FOR RENDERING LARGE COLUMNS ON THE HANDSONTABLE
function showFullString (event, title, longString) {
  $("#showLongStringDialog").html("<textarea class=\"code\" rows=\"10\" style=\"width: 100%\" readonly>"+ longString + "</textarea>")
  $("#showLongStringDialog").dialog({ title: title, position: {my: "left top", at: "left bottom", of: event.srcElement}})
}

function longColumnRenderer (instance, td, row, col, prop, value, cellProperties) {
  var longString = Handsontable.helper.stringify(value);
  

  if (longString.length>20) {
    //modify it
    longString = longString.substring(0, 20) + "<span>...<img class=\"font_size_images\" src=\"lib/resources/enlarge.png\" onclick=\"showFullString(event, '" + instance.getColHeader(col) + "', '"+longString+"')\"/></span>"
  }
  
  
  td.innerHTML = longString;
  return td;
}
//FUNCTIONS FOR RENDERING LARGE COLLUMNS ON THE HANDSONTABLE
//************************************************************






//***************************************************************
//MAIN FUNCTIONS
function buildPage(graphElements, DBGWAS_graph_tag2nodes)
{
    //this is basically main()
    $(function(){ // on dom ready
        $("#config button").prop("disabled", true); //disable the buttons

        //create the get fasta dialog
        dialogGetFASTANodes = $('#dialogGetFASTANodes').dialog({
            autoOpen: false,
            modal: true,
            maxWidth: 0.8*$(window).width(),
            maxHeight: 0.8*$(window).height()
        });

        //create the instructions dialog
        dialogInstructions = $('#dialogInstructions').dialog({
            autoOpen: false,
            modal: true,
            width: 0.8*$(window).width(),
            maxHeight: 0.8*$(window).height()
        });

        //create the cytoscape dialog
        cytoscapeExportDialog = $('#cytoscapeExportDialog').dialog({
            autoOpen: false,
            modal: true,
            width: 0.8*$(window).width(),
            maxHeight: 0.8*$(window).height()
        });

        //create cytoscape graph
        cy = cytoscape({

            container: document.getElementById('cy'),

            elements: graphElements,

            style: [ // the stylesheet for the graph
                {
                    selector: 'node',
                    style: {
                        'background-color': '#000000',
                        'label': 'data(info)',
                        'z-index': 10
                    }
                },

                {
                    selector: 'edge',
                    style: {
                        'width': 5,
                        'line-color': '#000000',
                        'target-arrow-color': '#000000',
                        'target-arrow-shape': 'triangle',
                        'font-size': '100',
                        'color': 'black',
                        'z-index': 10
                    }
                },
                {
                    selector: '.multiline-manual',
                    style: {
                        'text-wrap': 'wrap'
                    }
                },
                //style of the selection box
                {"selector":"core","style":{"selection-box-color":"#AAD8FF","selection-box-border-color":"#8BB0D0","selection-box-opacity":"0.5"}},
                //style of selected nodes
                {"selector":"node:selected","style":{"border-width":"20px","border-color":"#000000"}},
            ],
            boxSelectionEnabled: true,
            selectionType: 'additive'
        });

        //register what to do when user selects a node
        cy.nodes().on('select', function(evt){
            fillTable();
        });

        cy.nodes().on('unselect', function(evt){
            fillTable();
        });


        //add the helper for navigation
        cy.panzoom();

        //add the context Menus
        cy.contextMenus({
            menuItems: [
                {
                    id: 'select_all_nodes',
                    content: 'Select all nodes',
                    tooltipText: 'Select all nodes',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        selectAllNodes();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },


                {
                    id: 'unselect_all_nodes',
                    content: 'Unselect all nodes',
                    tooltipText: 'Unselect all nodes',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        unselectAllNodes();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },
                {
                    id: 'copy_fasta_to_clipboard',
                    content: 'Copy FASTA of selected nodes to clipboard',
                    tooltipText: 'Copy FASTA of selected nodes to clipboard',
                    selector: 'node, edge',
                    onClickFunction: function (event) {
                        copyFastaToClipboard();
                    },
                    hasTrailingDivider: true,
                    coreAsWell: true
                },
            ],
            menuItemClasses: ['custom-menu-item'],
            contextMenuClasses: ['custom-context-menu']

        });


        //populate the dropdown list for the tags
        Object.keys(DBGWAS_graph_tag2nodes).forEach(function(key) {
            $('#DBGWAS_graph_tag_Select').append($('<option>', {
                value: key,
                text: key,
            }));
        })

        //highlight the tag nodes
        $('#DBGWAS_graph_tag_Select').change(function(){
            selectNodesFromATag(this, cy, DBGWAS_graph_tag2nodes)
        });

        //say we are drawing the layout
        $("#PlWarning").html("Drawing layout...");

        //define the layout
        var layout = cy.layout({
            name: 'cytoscape-ngraph.forcelayout',
            async: {
// tell layout that we want to compute all at once:
                maxIterations: 2000,
                stepsPerCycle: 30,

// Run it till the end:
                waitForStep: false
            },
            physics: {
                /**
                 * Ideal length for links (springs in physical model).
                 */
                springLength: 300,

                /**
                 * Hook's law coefficient. 1 - solid spring.
                 */
                springCoeff: 0.00005,

                /**
                 * Coulomb's law coefficient. It's used to repel nodes thus should be negative
                 * if you make it positive nodes start attract each other :).
                 */
                gravity: -80,

                /**
                 * Theta coefficient from Barnes Hut simulation. Ranged between (0, 1).
                 * The closer it's to 1 the more nodes algorithm will have to go through.
                 * Setting it to one makes Barnes Hut simulation no different from
                 * brute-force forces calculation (each node is considered).
                 */
                theta: 1,

                /**
                 * Drag force coefficient. Used to slow down system, thus should be less than 1.
                 * The closer it is to 0 the less tight system will be.
                 */
                dragCoeff: 0,

                /**
                 * Default time step (dt) for forces integration
                 */
                timeStep: 20,
                iterations: 10000,
                fit: true,

                /**
                 * Maximum movement of the system which can be considered as stabilized
                 */
                stableThreshold: 0.000009
            },
            iterations: 10000,
            refreshInterval: 16, // in ms
            refreshIterations: 10, // iterations until thread sends an update
            stableThreshold: 2,
            animate: true,
            fit: true
        });

        //event when finishing the layout
        layout.on('layoutstop', function(){
            //say we are ready
            $("#PlWarning").html("Ready!");
            $("#config button").prop("disabled", false);
            window.callPhantom(); //tells phantom we are ready
        })

        //run the layout
        layout.run();

        //fills the instruction
        $("#dialogInstructions").html("\
    <ul>\
        <li>Navigation</li>\
        <ul>\
            <li>-Click and drag to move the screen;</li>\
            <li>-Click and drag a node to move it;</li>\
            <li>Use mouse wheel to zoom;</li>\
            <li>You can also use the navigation panel in the up left to navigate;</li>\
        </ul>\
        <li>Selecting a node</li>\
        <ul>      \
            <li>Press on a node to select it</li>\
            <li>Selecting a node will add it to the Node table in the bottom of the screen</li>\
            <li>To make a selection box, hold Ctrl and draw the box</li>\
            <li>You can easily select and unselect all nodes by right-clicking anywhere in the graph;</li>\
            <li>Press on a selected node to unselect it</li>\
            <li>Press anywhere in the graph to unselect all nodes</li>\
        </ul>\
    </ul>");

        //fills the cytoscape dialog instructions
        $("#cytoscapeExportDialog").html("\
          1. Download the graph <a download=\"graph2cytoscapeDesktop.json\" id=\"graph_cytoscape\"><b><u>here</u></b></a>, save it, and load it into Cytoscape (File > Import > Network > File...) <br/>\
          2. After the graph is loaded, download the style <a href=\"lib/xml/DBGWAS_cytoscape_style.xml\" download=\"DBGWAS_cytoscape_style.xml\"><b><u>here</u></b></a>, save it, and load it into Cytoscape (File > Import > Styles...) <br/>\
          3. To apply the style, go to the Style tab in the Control Panel and select DBGWAS_cytoscape_style.")

        //add the listeners to the download buttons in the cytoscape dialog instructions
        document.getElementById('graph_cytoscape').addEventListener('click', function () {
            var link = document.getElementById('graph_cytoscape');
            link.href = makeFile(JSON.stringify(cy.json()), cytoscapeDesktopGraph, 'application/json');
        }, false);


        //create the node table
        var tableSettings = {
            data: nodeTableData,
            columns: [
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                //TODO: NA temporarily removed
                //{type: 'text'},
                {type: 'text'},
                {renderer: longColumnRenderer},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {type: 'text'},
                {renderer: longColumnRenderer}
            ],
            colHeaders: [
                'Node ID',
                'Allele freq',
                'Pheno0',
                'Pheno1',
                //TODO: NA temporarily removed
                //'NA',
                'Annotation',
                'Significant?',
                'q-Value',
                'Est effect',
                'Wald stat',
                'Seq Length',
                'Sequence'
            ],
            copyColsLimit: 1000000,
            copyRowsLimitNumber: 1000000,
            readOnly: true,
            stretchH: 'all',
            wordWrap: false,
            allowInsertColumn: false,
            allowInsertRow: false,
            allowRemoveColumn: false,
            allowRemoveRow: false,
            autoColumnSize: {useHeaders: true},
            autoWrapCol: true,
            autoWrapRow: true,
            manualColumnResize: true,
            renderAllRows: true,
            width: widthTable,
            height: heightTable,
            columnSorting: true
        };

        var tableContainer = document.getElementById('nodeTable');
        var addedStyle = "width:" + widthTable.toString() + "px; " +
            "max-width:" + widthTable.toString() + "px; " +
            "height: " + heightTable.toString() + "px;"
        "max-height: " + heightTable.toString() + "px;"
        tableContainer.setAttribute("style", tableContainer.getAttribute("style") + addedStyle);
        table = new Handsontable(tableContainer, tableSettings);

        drawGradient();
        drawAlleles();
    }); // on dom ready
}
//MAIN FUNCTIONS
//***************************************************************