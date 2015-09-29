/* Handle web interface */



"use strict";



var ising = null;
var l = 32;
var tp = 2.3;

var timer_interval = 100; // in milliseconds
var ising_timer = null;
var mc_algorithm = "Metropolis";

var nstepspsmc = 1000; // number of steps per second for MC
var nstepspfmc = 100;  // number of steps per frame for MC
var mctot = 0.0;
var mcacc = 0.0;

var sum1 = 1e-300;
var sumU = 0;

var lnzref, eavref, cvref;



function getparams()
{
  l = get_int("L", 32);
  tp = get_float("temperature", 2.3);

  var arr = is2_exact(l, l, 1/tp);
  lnzref = arr[0];
  eavref = arr[1];
  cvref = arr[2];

  mc_algorithm = grab("mc_algorithm").value;

  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  mousescale = get_float("isingscale");
}



function changescale()
{
  mousescale = get_float("isingscale");
  paint();
}



function dometropolis()
{
  var istep, sinfo = "";
  var i, id, n = ising.n;

  ising.setproba(1.0/tp);
  //ising.em();
  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    for ( i = 0; i < n; i++ ) {
      mctot += 1.0;
      id = ising.pick();
      if ( ising.h <= 0 || rand01() < ising.proba[ising.h] ) {
        mcacc += 1;
        ising.flip(id);
      }
      sum1 += 1.0;
      sumU += ising.E;
    }
  }
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%, ";
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: '
         + roundto(sumU/sum1/n, 3) + " (Ref.: "
         + roundto(eavref/n, 3) + ").";
  return sinfo;
}



function dowolff()
{
  var istep, sinfo = "";
  var i, id, n = ising.n;
  var padd;

  padd = 1 - Math.exp(-2/tp);
  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    ising.wolff(padd);
    sum1 += 1.0;
    sumU += ising.E;
  }
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: '
         + roundto(sumU/sum1/n, 3) + " (Ref.: "
         + roundto(eavref/n, 3) + ").";
  return sinfo;
}



function paint()
{
  if ( !ising ) {
    return;
  }
  isingdraw(ising, "isingbox", mousescale);
}



function pulse()
{
  var sinfo;

  if ( mc_algorithm === "Metropolis" ) {
    sinfo = dometropolis();
  } else if ( mc_algorithm === "Wolff" ) {
    sinfo = dowolff();
  }
  grab("sinfo").innerHTML = sinfo;

  paint();
}



function stopsimul()
{
  if ( ising_timer !== null ) {
    clearInterval(ising_timer);
    ising_timer = null;
  }
  mctot = 0.0;
  mcacc = 0.0;
  sum1 = 1e-30;
  sumU = 0.0;
}



function pausesimul()
{
  if ( !ising ) {
    return;
  }
  if ( ising_timer !== null ) {
    clearInterval(ising_timer);
    ising_timer = null;
    grab("pause").value = "Resume";
  } else {
    ising_timer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").value = "Pause";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  ising = new Ising(l);
  installmouse("isingbox", "isingscale");
  ising_timer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



/* respond to critical parameter changes: restart simulation */
function changeparams()
{
  if ( ising_timer !== null ) {
    startsimul();
  }
}

