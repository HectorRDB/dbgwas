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