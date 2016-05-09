/* Handle web interface */



"use strict";



var ising = null;

var simtemp = null;
var ibeta = 0; // current temperature index
var tempfreq = 0.5; // tempering frequency
var ising_proba = null; // probability for the Metropolis algorithm
var ising_padd = null; // probability for the Wolff cluster algorithm
var hist = null; // energy histogram
var tpmctot = 1e-16, tpmcacc = 0;

var lnz_wham = null;

var timer_interval = 100; // in milliseconds
var ising_timer = null;
var mc_algorithm = "Metropolis";

var blocksizemc = 10; // number of steps for a block in the Metropolis algorithm
var nstepspsmc = 1000; // number of steps per second for MC
var nstepspfmc = 100;  // number of steps per frame for MC

var histplot = null;
var vplot = null;

function getparams()
{
  var l = get_int("L", 16);
  ising = new Ising(l);

  // set the temperature array
  var tpmin = get_float("tpmin", 1.5);
  var tpmax = get_float("tpmax", 3.1);
  var tpcnt = get_int("tpcnt", 17);
  var tp, dtp = (tpmax - tpmin) / tpcnt;
  var beta = newarr(tpcnt);
  var i;
  for ( i = 0; i < tpcnt; i++ ) {
    tp = tpmin + i * dtp;
    beta[i] = 1 / tp;
  }
  simtemp = new SimTemp(beta);
  ising_proba = newarr(tpcnt);
  ising_padd = newarr(tpcnt);
  hist = new Hist(tpcnt, -2*ising.n - 2, 2*ising.n + 1.999, 4);
  //console.log( hist.n, hist.xmin );

  for ( i = 0; i < tpcnt; i++ ) {
    var arr = is2_exact(l, l, beta[i]);
    simtemp.lnzref[i] = arr[0];
    simtemp.uref[i] = arr[1];
    simtemp.cvref[i] = arr[2];

    simtemp.lnz[i] = simtemp.lnzref[i];
    simtemp.hist[i] = 0;

    ising_proba[i] = newarr(5);
    var x = Math.exp(-4 * beta[i]);
    ising_proba[i][0] = 1.0;
    ising_proba[i][2] = x;
    ising_proba[i][4] = x * x;

    ising_padd[i] = 1 - Math.exp(-2*beta[i]);
  }
  ibeta = 0;

  lnz_wham = newarr(tpcnt);

  mc_algorithm = grab("mc_algorithm").value;

  blocksizemc = get_int("blocksizemc", 10);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  tpmcacc = 0;
  tpmctot = 1e-16;
  //mousescale = get_float("animationboxscale");
}



function changescale()
{
  //mousescale = get_float("animationboxscale");
  paint();
}



function dotempering()
{
  if ( rand01() < tempfreq ) {
    var jbeta = simtemp.jump(ibeta, ising.E, 1);
    tpmctot += 1;
    if ( ibeta != jbeta ) {
      tpmcacc += 1;
      ibeta = jbeta;
    }
  }
  simtemp.add(ibeta, ising.E);
  hist.arr[ibeta][(ising.E + ising.n*2)/4] += 1;
}


function ising_metropolis()
{
  var istep, j, id;

  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    for ( j = 0; j < blocksizemc; j++ ) {
      id = ising.pick();
      if ( ising.h <= 0
        || rand01() < ising_proba[ibeta][ising.h] ) {
        ising.flip(id);
      }
    }
    dotempering();
  }
}



function ising_wolff()
{
  var istep;

  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    ising.wolff(ising_padd[ibeta]);
    dotempering();
  }
}



function paint()
{
  if ( !ising ) {
    return;
  }
  // construct color from temperature
  var x = (ibeta + 0.5)/ simtemp.n;
  var color = rgb2str(Math.floor(204 * x), 0, Math.floor(204 * (1-x)));
  isingdraw(ising, "animationbox", 1.0, color);
}



function updatehistplot()
{
  var i, j, ntp = simtemp.n;
  var dat = "Energy";

  // prepare the header
  for ( j = 0; j < ntp; j++ ) {
    dat += ",Histogram " + (j+1);
  }
  dat += "\n";

  // refine the energy range
  var imin = 0;
  for ( ; imin <= ising.n; imin++ ) {
    for ( j = 0; j < ntp; j++ )
      if ( hist.arr[j][imin] > 0 ) break;
    if ( j < ntp ) break;
  }

  var imax = Math.floor(ising.n * 3 / 5);
  for ( ; imax >= imin; imax-- ) {
    for ( j = 0; j < ntp; j++ )
      if ( hist.arr[j][imax] > 0 ) break;
    if ( j < ntp ) break;
  }
  if ( imin > imax ) return;

  // fill in the energy histogram data
  for ( i = imin; i <= imax; i++ ) {
    var ep = -2 * ising.n + 4 * i;
    dat += "" + ep;
    for ( j = 0; j < ntp; j++ ) {
      dat += "," + hist.arr[j][i];
    }
    dat += "\n";
  }

  if ( histplot === null || 1 ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      xlabel: "<small>Energy</small>",
      ylabel: "<small>Histogram</small>",
      includeZero:true,
      axisLabelFontSize: 10,
      xRangePad: 2,
      plotter: barChartPlotter,
      width: w,
      height: h
    };
    histplot = new Dygraph(grab("histplot"), dat, options);
  } else {
    histplot.updateOptions({file: dat});
  }

  dat = null;
}



function updatevplot()
{
  var i, j, ntp = simtemp.n;
  // prepare the header
  var dat = "Temperature,lnZ<sub>WHAM</sub>,lnZ<sub>ref</sub>\n";

  // fill in the energy histogram data
  for ( j = 0; j < ntp; j++ ) {
    var tp = 1.0 / simtemp.beta[j];
    var lnzref = simtemp.lnzref[j] - simtemp.lnzref[0];
    dat += "" + tp + "," + lnz_wham[j] + ","
      + lnzref + "\n";
  }

  if ( vplot === null || 1 ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      xlabel: "<small>Temperature</small>",
      ylabel: "<small>ln Z</small>",
      includeZero:true,
      axisLabelFontSize: 10,
      xRangePad: 2,
      drawPoints: true,
      pointSize: 2,
      width: w,
      height: h
    };
    vplot = new Dygraph(grab("vplot"), dat, options);
  } else {
    vplot.updateOptions({file: dat});
  }
}



/* draw a color bar */
function drawcolorbar(target)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  var grd = ctx.createLinearGradient(0, 0, width, 0);
  grd.addColorStop(0, "#0000cc");
  grd.addColorStop(1, "#cc0000");

  ctx.fillStyle = grd;
  ctx.fillRect(0, 4, width, height-4);

  // draw grids
  if ( simtemp ) {
    var i, x;
    for ( i = 0; i < simtemp.n; i++ ) {
      x = width * i / (simtemp.n - 1.0);
      drawLine(ctx, x, 4, x, height, "#808080", 1);
    }
    x = width * ibeta / (simtemp.n - 1.0);
    drawLine(ctx, x, 0, x, height, "#000000", 3);
  }
}



function pulse()
{
  var sinfo;

  if ( mc_algorithm === "Metropolis" ) {
    ising_metropolis();
  } else if ( mc_algorithm === "Wolff" ) {
    ising_wolff();
  }
  drawcolorbar("tpscale");
  sinfo = "Temperature " + ibeta + "/" + simtemp.n + ", E " + ising.E
        + "<br>Tempering acc. ratio: " + (100.*tpmcacc/tpmctot).toPrecision(4) + "%";
  grab("sinfo").innerHTML = sinfo;

  paint();
  updatehistplot();
}



function stopsimul()
{
  if ( ising_timer !== null ) {
    clearInterval(ising_timer);
    ising_timer = null;
  }
  grab("pause").innerHTML = "&#9724;";
  tpmctot = 1e-16;
  tpmcacc = 0;
}



function pausesimul()
{
  if ( !ising ) {
    return;
  }
  if ( ising_timer !== null ) {
    clearInterval(ising_timer);
    ising_timer = null;
    grab("pause").innerHTML = "&#10704";
  } else {
    ising_timer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").innerHTML = "&#9724;";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  //installmouse("animationbox", "animationboxscale");
  ising_timer = setInterval(
    function(){ pulse(); },
    timer_interval );
}



function pausesimul2()
{
  // skip a mouse-move
  //if ( mousemoved > 0 ) {
  //  return;
  //}
  if ( !ising ) {
    startsimul();
  //} else if ( mousemoved === 0 ) {
  } else {
    pausesimul();
  }
}


function runwham()
{
  if ( ising_timer ) {
    pausesimul();
  }
  if ( simtemp ) {
    var wham_method = grab("WHAM-method").value;
    var wham_itmax = get_int("WHAM-itmax", 1000);
    var wham_nlogtol = get_int("WHAM-tol", -7);
    var wham_tol = Math.pow(10, wham_nlogtol);
    var wham_damp = get_float("WHAM-damp", 1.0);
    var mdiis_nbases = get_int("MDIIS-nbases", 5);
    var mdiis_method = grab("MDIIS-method").value;
    var kth_threshold = get_float("KTH-threshold", 10.0);

    var err = whamx(hist, simtemp.beta, lnz_wham, 0, null,
        wham_damp, mdiis_nbases, mdiis_method,
        kth_threshold, 0, wham_itmax,
        wham_tol, 0, wham_method);
    grab("sinfo").innerHTML = err;
    updatevplot();
    showtab("wham-params");
  }
}


/* respond to critical parameter changes: restart simulation */
function changeparams()
{
  if ( ising_timer !== null ) {
    startsimul();
  }
}



function showtab(who)
{
  who = grab(who);
  var par = who.parentNode;
  var c = par.childNodes;
  var i, iwho, k = 0;

  // arrange the tabs
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      if ( c[i] !== who ) {
        c[i].style.zIndex = k;
        k += 1;
      } else {
        iwho = k;
      }
    }
  }
  who.style.zIndex = k;

  // arrange the clickable tab titles
  k += 1;
  var pt = grab("tabsrow");
  pt.style.zIndex = k;
  var ct = pt.childNodes, ik = 0;
  for ( i = 0; i < ct.length; i++ ) {
    if ( ct[i].tagName ) {
      if ( ik === iwho ) {
        ct[i].style.fontWeight = "bold";
        ct[i].style.borderTop = "2px solid #c0c0d0";
      } else {
        ct[i].style.fontWeight = "normal";
        ct[i].style.borderTop = "0px solid #e0e0f0";
      }
      ik += 1;
    }
  }
}



function resizecontainer(a)
{
  var canvas = grab("animationbox");
  var ctx = canvas.getContext("2d");
  var w, h;
  if ( a === null || a === undefined ) {
    w = canvas.width;
    h = canvas.height;
  } else {
    a = parseInt( grab(a).value );
    w = h = a;
    canvas.width = w;
    canvas.height = h;
  }
  ctx.font = "24px Verdana";
  ctx.fillText("Click to start", w/2-40, h/2-10);

  var hsbar = 30; // height of the global scaling bar
  var hcbar = 40; // height of the control bar
  var htbar = 30; // height of the tabs bar
  var wr = h*3/4; // width of the plots
  var wtab = 560; // width of the tabs
  var htab = 280;

  grab("simulbox").style.width = "" + w + "px";
  grab("simulbox").style.height = "" + h + "px";
  grab("simulbox").style.top = "" + hsbar + "px";
  grab("controlbox").style.top = "" + (h + hsbar) + "px";
  //grab("animationboxscale").style.width = "" + (w - 100) + "px";
  grab("tpscale").style.width = "" + (w - 220) + "px";
  drawcolorbar("tpscale");
  histplot = null;
  grab("histplot").style.left = "" + w + "px";
  grab("histplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + hcbar + "px";
  grab("histplot").style.height = "" + h/2 + "px";
  vplot = null;
  grab("vplot").style.left = "" + w + "px";
  grab("vplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + (h/2 + hcbar) + "px";
  grab("vplot").style.height = "" + h/2 + "px";
  grab("tabsrow").style.top = "" + (h + hsbar + hcbar) + "px";
  grab("tabsrow").style.width = "" + wtab + "px";

  var c = grab("container").childNodes;
  var i;
  // tabs
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      c[i].style.top = "" + (h + hsbar + hcbar + htbar) + "px";
      c[i].style.width = "" + (wtab - 20) + "px";
      c[i].style.height = "" + htab + "px";
    }
  }
  grab("sinfo").style.top = "" + (h + hsbar + hcbar + htbar) + "px";
  grab("sinfo").style.left = "" + (wtab + 10) + "px";
  grab("sinfo").style.width = "" + (w + wr - wtab - 20) + "px";
  grab("sinfo").innerHTML = "simulation information";

  grab("container").style.height = "" + (h + hsbar + hcbar + htbar + htab) + "px";
  grab("container").style.width = "" + (w + wr) + "px";
}


function init()
{
  resizecontainer();
  showtab("system-params");
}
