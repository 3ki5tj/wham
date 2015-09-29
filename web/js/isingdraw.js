


"use strict";



/* draw all atoms in the box
 * `target` is the canvas  */
function isingdraw(ising, target, userscale)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + 1
  var dx = 1.0 * width / (ising.l + 1) * userscale;
  var dy = 1.0 * height / (ising.l + 1) * userscale;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // draw each spin
  var l = ising.l, id = 0;
  for ( var i = 0; i < l; i++ ) {
    for ( var j = 0; j < l; j++, id++ ) {
      var x = (i + 0.5 - 0.5*l) * dx + width * 0.5;
      var y = (j + 0.5 - 0.5*l) * dy + height * 0.5;
      var radius = dx * 0.5;
      var spotcolor, color;
      if ( ising.s[id] > 0 ) {
        color = "#202020";
        spotcolor = "#808080";
      } else {
        color = "#aaaaaa";
        spotcolor = "#ffffff";
      }
      paintBall(ctx, x, y, radius, color, spotcolor);
    }
  }
}



