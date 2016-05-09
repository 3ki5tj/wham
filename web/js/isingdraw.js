


"use strict";



/* draw all atoms in the box
 * `target` is the canvas  */
function isingdraw(ising, target, userscale,
    ballColor)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + 1
  var dx = 1.0 * width / (ising.l) * userscale;
  var dy = 1.0 * height / (ising.l) * userscale;
  var dx1 = Math.floor(dx - 1), dy1 = Math.floor(dy - 1);
  if ( dx1 < 1 ) dx1 = 1;
  if ( dy1 < 1 ) dy1 = 1;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // draw each spin
  var l = ising.l, id = 0;
  var upColor, downColor;
  if ( !ballColor ) {
    upColor = "#202020";
    downColor = "#cccccc";
  } else {
    upColor = ballColor;
    downColor = "#cccccc";
  }

  var grd_up = upColor;
  var grd_down = downColor;

  for ( var i = 0; i < l; i++ ) {
    for ( var j = 0; j < l; j++, id++ ) {
      var x = Math.floor( (i - 0.5*l) * dx + width * 0.5 );
      var y = Math.floor( (j - 0.5*l) * dy + height * 0.5 );
      var spotcolor, color;
      if ( ising.s[id] > 0 ) {
        color = upColor;
        ctx.fillStyle = grd_up;
      } else {
        color = downColor;
        ctx.fillStyle = grd_down;
      }
      ctx.fillRect(x, y, dx1, dy1)
    }
  }
}



