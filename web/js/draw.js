/* graphics routines */



"use strict";



/* draw line */
function drawLine(ctx, xi, yi, xj, yj, color, lineWidth)
{
  if ( color  ) {
    ctx.strokeStyle = color;
  }
  if ( lineWidth ) {
    ctx.lineWidth = lineWidth;
  }
  ctx.beginPath();
  ctx.moveTo(xi, yi);
  ctx.lineTo(xj, yj);
  ctx.stroke();
}



// draw a 3-layer line
function drawLineGradient(ctx, xi, yi, xj, yj, grd)
{
  if ( !grd ) {
    grd = [ { color: '#aaaaaa', width:8 },
            { color: '#bbbbbb', width:4 },
            { color: '#cccccc', width:2 } ];
  }

  for ( var i = 0; i < grd.length; i++ ) {
    drawLine(ctx, xi, yi, xj, yj, grd[i].color, grd[i].width);
  }
}



/* draw a ball that is centered at (x, y) with radius r
 * color is the color of the ball
 * the format of color is "#rrggbb" */
function drawBall(ctx, x, y, r, color, lineWidth)
{
  if ( color  ) {
    ctx.strokeStyle = color;
  }
  if ( lineWidth ) {
    ctx.lineWidth = lineWidth;
  }
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI);
  ctx.closePath();
  ctx.stroke();
}



/* draw a ball that is centered at (x, y) with radius r
 * color is the color of the ball
 * spotcolor is the color of the spotlight
 * the format of color is "#rrggbb" */
function paintBall(ctx, x, y, r, color, spotcolor,
    spotx, spoty, spotr)
{
  if ( spotcolor === undefined || spotcolor === null ) {
    spotcolor = "#a0a0a0";
  }
  if ( spotx === undefined || spotcolor === null ) {
    spotx = r * 0.3;
  }
  if ( spoty === undefined || spoty === null ) {
    spoty = r * 0.4;
  }
  if ( spotr === undefined || spotr === null ) {
    spotr = r * 0.1;
  }
  var grd = ctx.createRadialGradient(x + spotx, y - spoty, spotr, x, y, r);
  grd.addColorStop(0, spotcolor); // spotlight color
  grd.addColorStop(1, color); // ball color
  ctx.fillStyle = grd;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI);
  ctx.closePath();
  ctx.fill();
}



/* compute the contact point to start a bond from i to j */
function getContactPoint(xi, xj, radius)
{
  var rji, xji = newarr(D);

  rji = vdistx(xji, xj, xi);
  vsmul(xji,  radius / rji);
  return vinc(xji, xi);
}



/* convert RGB values to string */
function rgb2str(r, g, b)
{
  r = Math.floor(r).toString(16);
  if ( r.length == 1 ) {
    r = "0" + r;
  }

  g = Math.floor(g).toString(16);
  if ( g.length == 1 ) {
    g = "0" + g;
  }

  b = Math.floor(b).toString(16);
  if ( b.length == 1 ) {
    b = "0" + b;
  }

  return "#" + r + g + b;
}



/* get the red, green and blue components of a color string */
function parseRGB(color)
{
  if ( color == undefined ) {
    console.log("");
  }
  if ( color.substr(0, 3) === "rgb" ) {
    var i0 = color.indexOf("(") + 1;
    var i1 = color.lastIndexOf(")");
    var s = color.substring(i0, i1).split(",");
    return {
      r: parseInt(s[0], 10),
      g: parseInt(s[1], 10),
      b: parseInt(s[2], 10)
    };
  } else {
    return {
      r: parseInt(color.substr(1, 2), 16),
      g: parseInt(color.substr(3, 2), 16),
      b: parseInt(color.substr(5, 2), 16)
    };
  }
}



/* add transparency to color */
function transpColor(color, transp)
{
  var c = parseRGB(color);
  return "rgba(" + c.r + ", " + c.g + ", " + c.b + ", " + transp + ")";
}



function randHueColor(cmin, cmax)
{
  var x = Math.random() * 6;
  var i = Math.floor( x ), r = 0, g = 0, b = 0;

  if ( cmin === undefined || cmin === null ) {
    cmin = 0;
  }
  if ( cmax === undefined || cmax === null ) {
    cmax = 255;
  }
  var cvar = cmax - cmin + 1;
  x -= i;
  if ( i < 1 ) { // red to yellow
    r = cmax;
    g = cmin + Math.floor( cvar * x );
  } else if ( i < 2 ) { // yellow to green
    r = cmin + Math.floor( cvar * (1 - x) );
    g = cmax;
  } else if ( i < 3 ) { // green to cyan
    g = cmax;
    b = cmin + Math.floor( cvar * x );
  } else if ( i < 4 ) { // cyan to blue
    g = cmin + Math.floor( cvar * (1 - x) );
    b = cmax;
  } else if ( i < 5 ) { // blue to magenta
    b = cmax;
    r = cmin + Math.floor( cvar * x );
  } else {
    b = cmin + Math.floor( cvar * (1 - x) );
    r = cmax;
  }
  return rgb2str(r, g, b);
}



/* darken a color
 * fac === 1.0 means to preserve the original order */
function darkenColor(colorStr, fac)
{
  if ( !fac ) {
    fac = 0.5;
  }
  var color = parseRGB(colorStr);
  color.r = Math.floor(color.r * fac);
  color.g = Math.floor(color.g * fac);
  color.b = Math.floor(color.b * fac);
  return rgb2str(color.r, color.g, color.b);
}



/* lighten a color
 * fac === 1.0 means to preserve the original order */
function lightenColor(colorStr, fac)
{
  if ( !fac ) {
    fac = 0.5;
  }
  var color = parseRGB(colorStr);
  color.r = Math.floor(color.r * fac + (1 - fac) * 255);
  color.g = Math.floor(color.g * fac + (1 - fac) * 255);
  color.b = Math.floor(color.b * fac + (1 - fac) * 255);
  return rgb2str(color.r, color.g, color.b);
}



/* get the greyscale version of a color */
function greyColor(colorStr, fac)
{
  var color = parseRGB(colorStr);
  var x = Math.floor((color.r + color.g + color.b)/3);
  return rgb2str(x, x, x);
}



/* bar chart for Dygraph */
function barChartPlotter(e)
{
  var ctx = e.drawingContext;
  var points = e.points;
  var y_bottom = e.dygraph.toDomYCoord(0);
  var i;

  ctx.fillStyle = transpColor( e.color, 0.5 );

  // Find the minimum separation between x-values.
  // This determines the bar width.
  var min_sep = Infinity;
  for ( i = 1; i < points.length; i++ ) {
    var sep = points[i].canvasx - points[i - 1].canvasx;
    if (sep < min_sep) {
      min_sep = sep;
    }
  }
  var bar_width = Math.floor(2.0 / 3 * min_sep);

  // Do the actual plotting.
  for ( i = 0; i < points.length; i++ ) {
    var p = points[i];
    var center_x = p.canvasx;

    ctx.fillRect(center_x - bar_width / 2, p.canvasy,
        bar_width, y_bottom - p.canvasy);

    ctx.strokeRect(center_x - bar_width / 2, p.canvasy,
        bar_width, y_bottom - p.canvasy);
  }
}
