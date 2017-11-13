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

//Render URL to file

var system = require("system");
/*
Render given url
*/
RenderUrlToFile = function(url, outFile) {
  webpage = require("webpage");
  page = webpage.create();
  page.viewportSize = {
      width: 2560,
      height: 1600
  };
  page.onError = function() {
    console.log("Error on the webpage!");
    phantom.exit(1);
  }
  page.onCallback = function() {
      //get the bounding box of the cy element
      var bb = page.evaluate(function () { 
          return document.getElementById('cy').getBoundingClientRect(); 
      });

      //clip the page to this bb
      page.clipRect = {
          top:    bb.top,
          left:   bb.left,
          width:  bb.width,
          height: bb.height
      };

      //render it
      page.render(outFile);
      page.close();
      phantom.exit();
  };

  page.open("file:///" + url, function(status) {
    if (status !== "success") {
      console.log("Unable to render '" + url + "'");
      phantom.exit(1);
    }
  });
}

RenderUrlToFile(system.args[1], system.args[2])