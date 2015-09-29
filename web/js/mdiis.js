/* modified direct inversion of the iterative subspace (MDIIS) method */




"use strict";





/* open an mdiis object */
function MDIIS(npt, mnb, getres, obj, verbose)
{
  this.npt = npt;
  this.mnb = mnb; // maximal size of the basis
  this.nb = 0;
  this.ibQ = -1;
  this.getres = getres; // function to get residues
  this.obj = obj; // handle to the object
  var mnb1 = this.mnb + 1;
  this.f = newarr2d(mnb1, npt); // basis
  this.res = newarr2d(mnb1, npt); // residual vectors
  this.mat = newarr(mnb * mnb); // correlations of residues
  this.mat2 = newarr(mnb1 * mnb1); // temporary matrix for LU decomposition
  this.coef = newarr(mnb1); // coefficients of combining base vectors
  this.fbest = newarr(npt);
  this.errmin = 1e300;
  this.verbose = verbose;
  this.TINY = 1e-20;
}



/* solve the coefficients of combination */
MDIIS.prototype.solve = function()
{
  var nb = this.nb, nb1 = this.nb + 1, mnb = this.mnb, i, j;

  // right-hand side of the equation
  for ( i = 0; i < nb; i++ ) {
    this.coef[i] = 0;
  }
  this.coef[nb] = -1;

  // copy the matrix, for the content is to be destroyed
  for ( i = 0; i < nb; i++ ) {
    for ( j = 0; j < nb; j++ ) {
      this.mat2[i*nb1 + j] = this.mat[i*mnb + j];
    }
  }
  for ( i = 0; i < nb1; i++ ) {
    this.mat2[i*nb1 + nb] = -1;
    this.mat2[nb*nb1 + i] = -1;
  }
  this.mat2[nb*nb1 + nb] = 0;

  if ( lusolve(this.mat2, this.coef, nb1, this.TINY) != 0 ) {
    console.log("MDIIS lusolve failed");
    return -1;
  }
  return 0;
}



/* construct the new f */
MDIIS.prototype.gen = function(f, normalize, damp)
{
  var ib, il, npt = this.npt, nb = this.nb;

  if ( damp === undefined ) {
    damp = 1;
  }

  for ( il = 0; il < npt; il++ )
    this.f[nb][il] = 0;
  for ( ib = 0; ib < nb; ib++ ) {
    var coef = this.coef[ib];
    for ( il = 0; il < npt; il++ ) {
      this.f[nb][il] += coef * (this.f[ib][il] + damp * this.res[ib][il]);
    }
  }

  if ( normalize ) {
    normalize(this.f[nb], npt);
  }

  // f = this.f[nb]
  cparr(f, this.f[nb], npt);
}



/* compute the dot product */
function mdiis_getdot(a, b, n)
{
  var i, x = 0;

  for ( i = 0; i < n; i++ ) {
    x += a[i] * b[i];
  }
  return x / n;
}



/* build the residue correlation matrix */
MDIIS.proto.build = function(f, res)
{
  var i, npt = this.npt;

  this.nb = 1;

  for ( i = 0; i < npt; i++ ) {
    this.f[0][i] = f[i];
    this.res[0][i] = res[i];
  }

  this.mat[0] = mdiis_getdot(this.res[0], this.res[0], npt);

  this.ibQ = -1;
  return 0;
}



/* try to add the new vector `f` and its residue `res`
 * into the base using the Kovalenko-Ten-no-Hirata scheme */
MDIIS.prototype.update_kth = function(f, res, err, threshold)
{
  var i, ibmin, ib, npt = this.npt;
  var dot, min;

  var nb = this.nb;
  var mnb = this.mnb;

  // save this function if it achieves the minimal error so far
  if ( err < this.errmin ) {
    cparr(this.fbest, this.f[nb], npt);
    this.errmin = err;
  }

  // choose the base with the smallest residue
  ibmin = 0;
  for ( i = 1; i < nb; i++ ) {
    // the diagonal represents the error
    if ( this.mat[i*mnb + i] < this.mat[ibmin*mnb + ibmin] ) {
      ibmin = i;
    }
  }
  min = this.mat[ibmin*mnb + ibmin];
  dot = mdiis_getdot(res, res, npt);

  if ( dot > threshold * threshold * min ) {
    // rebuild the basis
    this.build(this.f[ibmin], this.res[ibmin]);
    return 0;
  }

  if ( nb < mnb ) {
    ib = nb;
    this.nb = ++nb;
  } else {
    ib = (this.ibQ + 1) % mnb;
    this.ibQ = ib;
  }

  // replace base ib by f
  for ( i = 0; i < npt; i++ ) {
    this.f[ib][i] = f[i];
    this.res[ib][i] = res[i];
  }

  // update the residue correlation matrix
  // note: we do not need to update the last row & column
  for ( i = 0; i < nb; i++ ) {
    this.mat[i*mnb + ib] = this.mat[ib*mnb + i]
      = mdiis_getdot(this.res[i], res, npt);
  }
  return ib;
}



/* try to add the new vector `f` and its residue `res`
 * into the base using the Howard-Pettitt scheme */
MDIIS.prototype.update_hp = function(f, res, err)
{
  var i, ibmax, ib, nb = this.nb, mnb = this.mnb, npt = this.npt;

  // save this function if it achieves the minimal error so far
  if ( err < this.errmin ) {
    cparr(this.fbest, this.f[nb], npt);
    this.errmin = err;
  }

  // find the base with the largest residue
  ibmax = 0;
  for ( i = 1; i < nb; i++ ) {
    // the diagonal represents the error
    if ( this.mat[i*mnb + i] > this.mat[ibmax*mnb + ibmax] ) {
      ibmax = i;
    }
  }

  // if we are updating the same vector from the previous step
  // then we are stuck, so rebuild the basis
  // This condition cannot be true until we have a full basis
  if ( ibmax == this.ibQ ) {
    // rebuild the basis from f
    // if we rebuild the basis from f[ibmin]
    // it is more likely to enter a limit cycle
    this.build(f, res);
    return 0;
  }

  if ( nb < this.mnb ) {
    ib = nb;
    this.nb = ++nb;
  } else {
    ib = ibmax;
    // Note: we only set ibQ if the basis is full
    this.ibQ = ib;
  }

  // replace base ib by f
  for ( i = 0; i < npt; i++ ) {
    this.f[ib][i] = f[i];
    this.res[ib][i] = res[i];
  }

  // update the residue correlation matrix
  // note: we do not need to update the last row & column
  for ( i = 0; i < nb; i++ ) {
    this.mat[i*mnb + ib] = this.mat[ib*mnb + i]
      = mdiis_getdot(this.res[i], res, npt);
  }
  return ib;
}



/* try to add the new vector `f` and its residue `res`
 * into the base using the Howard-Pettitt scheme */
MDIIS.prototype.update_hpl = function(f, res, err)(
{
  var i, ibmax, ib, nb = this.nb, mnb = this.mnb, npt = this.npt;

  // save this function if it achieves the minimal error so far
  if ( err < this.errmin ) {
    cparr(this.fbest, this.f[nb], npt);
    this.errmin = err;
  }

  // find the base with the largest (except the previous) residue
  ibmax = ( this.ibQ == 0 && nb > 1 ) ? 1 : 0;
  for ( i = ibmax + 1; i < nb; i++ ) {
    if ( i == this.ibQ ) continue;
    // the diagonal represents the error
    if ( this.mat[i*mnb + i] > this.mat[ibmax*mnb + ibmax] ) {
      ibmax = i;
    }
  }

  if ( nb < mnb ) {
    ib = nb;
    this.nb = ++nb;
  } else {
    ib = ibmax;
    // Note: we only set ibQ if the basis is full
    this.ibQ = ib;
  }

  // replace base ib by f
  for ( i = 0; i < npt; i++ ) {
    this.f[ib][i] = f[i];
    this.res[ib][i] = res[i];
  }

  // update the residue correlation matrix
  // note: we do not need to update the last row & column
  for ( i = 0; i < nb; i++ ) {
    this.mat[i*mnb + ib] = this.mat[ib*mnb + i]
      = mdiis_getdot(this.res[i], res, npt);
  }
  return ib;
}



/* try to add the new vector `f` and its residue `res`
 * into the base */
MDIIS.prototype.update = function(f, res, err)
{
  var i, ib, jb, nb = this.nb, mnb = this.mnb, npt = this.npt;
  var dot, max;

  // save this function if it achieves the minimal error so far
  if ( err < this.errmin ) {
    cparr(this.fbest, this.f[nb], npt);
    this.errmin = err;
  }

  // choose the base with the largest residue
  ib = 0;
  for ( i = 1; i < nb; i++ ) {
    // the diagonal represents the error
    if ( this.mat[i*mnb + i] > this.mat[ib*mnb + ib] ) {
      ib = i;
    }
  }
  max = this.mat[ib*mnb + ib];

  dot = mdiis_getdot(res, res, npt);

  if ( dot >= max ) {
    if ( nb >= 2 ) {
      // shrink the basis by removing the base with
      // the largest residue and try again
      jb = nb - 1;
      if ( ib != jb ) {
        // move the last base to position ib
        for ( i = 0; i < npt; i++ ) {
          this.f[ib][i] = this.f[jb][i];
          this.res[ib][i] = this.res[jb][i];
        }

        for ( i = 0; i < nb - 1; i++ ) {
          if ( i == ib ) continue;
          this.mat[i*mnb + ib] = this.mat[i*mnb + jb];
          this.mat[ib*mnb + i] = this.mat[jb*mnb + i];
        }
        this.mat[ib*mnb + ib] = this.mat[jb*mnb + jb];
      }
      this.nb--;
      return -1;
    } else {
      /* rebuild the basis from f */
      this.build(f, res);
      return 0;
    }
  }

  if ( nb < mnb ) {
    ib = nb;
    this.nb = ++nb;
  }

  /* replace base ib by f */
  for ( i = 0; i < npt; i++ ) {
    this.f[ib][i] = f[i];
    this.res[ib][i] = res[i];
  }

  /* update the residue correlation matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ ) {
    this.mat[i*mnb + ib] = this.mat[ib*mnb + i]
      = mdiis_getdot(this.res[i], res, npt);
  }
  return ib;
}



function iter_mdiis(f, npt, getres, normalize, obj,
    nbases, damp, update_method, threshold,
    itmin, itmax, tol, verbose)
{
  var it, ibp = 0, ib, success;
  var err, errp;

  // open an mdiis object
  var mdiis = new MDIIS(npt, nbases, getres, obj, verbose);
  // use the space of the last array for the current residue
  var res = mdiis.res[mdiis.mnb];

  // construct the initial base set
  mdiis.errmin = err = errp = mdiis.getres(obj, f, res);
  mdiis.build(f, res);

  for ( it = 0; it < itmax; it++ ) {
    // obtain a set of optimal coefficients of combination
    mdiis.solve();
    /* generate a new f from the set of coefficients */
    mdiis.gen(f, normalize, damp);
    err = mdiis.getres(obj, f, res);
    /* add the new f into the basis */
    if ( update_method == "KTH" ) {
      ib = mdiis.update_kth(f, res, err, threshold);
    } else if ( update_method == "HP" ) {
      ib = mdiis.update_hp(f, res, err);
    } else if ( update_method == "HPL" ) {
      ib = mdiis.update_hpl(f, res, err);
    } else {
      ib = mdiis.update(f, res, err);
    }

    if ( ib >= 0 ) {
      ibp = ib;
      errp = err;
    }

    if ( err < tol && it >= itmin ) {
      break;
    }
  }

  if ( err < tol ) {
    success = 1;
  } else { /* use the backup version */
    success = 0;
    cparr(f, mdiis.fbest, npt);
    err = mdiis.getres(obj, f, res);
  }

  console.log("MDIIS finished in ", it + 1, " steps, error ", err,
      " (", (success ? "succeeded" : "failed"), ")");
  return err;
}



