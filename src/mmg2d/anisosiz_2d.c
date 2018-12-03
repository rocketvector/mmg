/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/
/**
 * \file mmg2d/anisosiz_2d.c
 * \brief Interpolation of metrics
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param k elt index
 * \param i vertex index in triangle \a k
 *
 * \return 0 if fail, 1 if success
 *
 * Impose default metric (isotropic, with size hmax) at vertex i in triangle k
 * (don't take into account the local parameters). Set the point flag to 1 to be
 * able to truncate it with the local params later.
 *
 */
int MMG2D_defaultmet_2d(MMG5_pMesh mesh,MMG5_pSol met,int k,char i) {
  MMG5_pTria       pt;
  double           *m,isqhmax;
  int              ip;

  isqhmax = mesh->info.hmax;

  isqhmax = 1.0 / (isqhmax*isqhmax);
  pt = &mesh->tria[k];
  ip = pt->v[i];
  m = &met->m[3*ip];

  m[0] = isqhmax;
  m[1] = 0.0;
  m[2] = isqhmax;

  mesh->point[ip].flag = 1;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param k index of the tria in which we work
 * \param i index of the point on which we want to compute the metric
 *
 * \return 1 if success, 0 if fail
 *
 * Calculate anisotropic metric tensor at (boundary) vertex i in triangle k on
 * account of geometric approximation of the corresponding curve (taking into
 * account the local parameters). Set the point flag to 2 to ignore it whem
 * imposing the local parameters later.
 *
 */
int MMG2D_defmetbdy_2d(MMG5_pMesh mesh,MMG5_pSol met,int k,char i) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0,p1,p2;
  MMG5_pPar       ppa;
  double          hausd,hmin,hmax,sqhmin,sqhmax,ux,uy,ll,li,ps1,ps2,lm,ltmp,pv;
  double          M1,M2,t1[2],t2[2],b1[2],b2[2],*n,*m;
  double          gpp1[2],gpp2[2];
  int             ilist,iel,ip,ip1,ip2,it[2],l,list[MMG2D_LONMAX+2];
  int8_t          isloc,hausdloc;
  char            i0,i1,i2,j;
  static int8_t   mmgWarn0=0,mmgWarn1=0,mmgWarn2=0;

  hmin   = mesh->info.hmin;
  hmax   = mesh->info.hmax;
  hausd  = mesh->info.hausd;

  pt = &mesh->tria[k];
  ip = pt->v[i];
  p0 = &mesh->point[ip];
  m = &met->m[3*ip];

  ip1 = ip2 = 0;
  ilist = MMG2D_boulet(mesh,k,i,list);

  /* Local parameters if needed: note that the hausdorff param is only looked if
   * imposed on an edge */
  isloc = 0;
  hausdloc = 0;
  if ( mesh->info.npar ) {
    /* Minimum size feature imposed by triangles */
    for (k=0; k<ilist; k++) {
      iel = list[k]/3;
      pt = &mesh->tria[iel];
      assert ( MG_EOK(pt) );

      /* Retrieve local parameters associated to triangle k */
      for (l=0; l<mesh->info.npar; l++) {
        ppa = &mesh->info.par[l];
        if ( ppa->elt == MMG5_Triangle && ppa->ref == pt->ref ) {
          if ( !isloc ) {
            hmin = ppa->hmin;
            hmax = ppa->hmax;
            isloc = 1;
          }
          else {
            hmin  = MG_MAX ( hmin, ppa->hmin );
            hmax  = MG_MIN ( hmax, ppa->hmax );
          }
          break;
        }
      }
      /* Minimum size feature imposed by the boundary edge */
      for ( i=0; i<3; i++ ) {
        if ( !MG_EDG(pt->tag[i]) ) continue;

        if ( mesh->info.npar ) {
          for (l=0; l<mesh->info.npar; l++) {
            ppa = &mesh->info.par[l];
            if ( ppa->elt == MMG5_Edg && ppa->ref == pt->edg[i] ) {
              if ( !hausdloc ) {
                hausd = ppa->hausd;
                hausdloc = 1;
              }
              else {
                hausd = MG_MIN(hausd,ppa->hausd);
              }
              if ( !isloc ) {
                hmax = ppa->hmax;
                hmin = ppa->hmin;
              }
              else {
                hmin  = MG_MAX ( hmin, ppa->hmin );
                hmax  = MG_MIN ( hmax, ppa->hmax );
              }
              break;
            }
          }
        }
      }
    }

    /* Minimum size feature imposed by the vertex */
    for (l=0; l<mesh->info.npar; l++) {
      ppa = &mesh->info.par[l];
      if ( ppa->elt == MMG5_Vertex && ppa->ref == p0->ref ) {
        if ( !isloc ) {
          hmin = ppa->hmin;
          hmax = ppa->hmax;
          isloc = 1;
        }
        else {
          hmin  = MG_MAX ( hmin, ppa->hmin );
          hmax  = MG_MIN ( hmax, ppa->hmax );
        }
        break;
      }
    }
  }

  if ( hmin > hmax ) {
    if ( !mmgWarn2 ) {
      assert ( isloc && "Non compatible local parameters" );
      fprintf(stderr,"\n  ## Warning: %s: Non compatible local parameters:\n"
              " hmin (%.15lg) > hmax (%.15lg).\nhmax ignored.",__func__,hmin,hmax);
      hmax = MMG5_HMINMAXGAP*hmin;
    }
    mmgWarn2 = 1;
  }
  sqhmin   = hmin*hmin;
  sqhmax   = hmax*hmax;

  /* Recover the two boundary edges meeting at ip */
  for (l=0; l<ilist; l++) {
    iel = list[l] / 3;
    pt = &mesh->tria[iel];

    i0 = list[l] % 3;
    i1 = MMG5_inxt2[i0];
    i2 = MMG5_iprv2[i0];

    if ( MG_EDG(pt->tag[i1]) ) {
      if ( ip1 == 0 ) {
        ip1 = pt->v[i2];
        it[0] = 3*iel+i1;
      }
      else if ( ip1 != pt->v[i2] ) {
        if ( ip2 == 0 ) {
          ip2 = pt->v[i2];
          it[1] = 3*iel+i1;
        }
        else if ( ip2 != pt->v[i2] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the"
                    " intersection of 3 edges. abort.\n",__func__);
          }
          return 0;
        }
      }
    }

    if ( MG_EDG(pt->tag[i2]) ) {
      if ( ip1 == 0 ) {
        ip1 = pt->v[i1];
        it[0] = 3*iel+i2;
      }
      else if ( ip1 != pt->v[i1] ) {
        if ( ip2 == 0 ) {
          ip2 = pt->v[i1];
          it[1] = 3*iel+i2;
        }
        else if ( ip2 != pt->v[i1] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the"
                    " intersection of 3 edges. abort.\n",__func__);
          }
          return 0;
        }
      }
    }
  }

  /* Check that there are exactly two boundary points connected at p0 */
  if ( ip1 == 0 || ip2 == 0 ) {
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 point that is not"
              "at the intersection of 2 edges. abort.\n",__func__);
    }
    return 0;
  }

  lm = sqhmax;

  /* Curvature of both boundary edges meeting at ip */
  for (j=0; j<2; j++) {
    iel = it[j] / 3;
    pt = &mesh->tria[iel];
    i0 = it[j] % 3;
    i1 = MMG5_inxt2[i0];
    i2 = MMG5_iprv2[i0];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];

    p1 = &mesh->point[ip1];
    p2 = &mesh->point[ip2];

    ux = p2->c[0] - p1->c[0];
    uy = p2->c[1] - p1->c[1];
    ll = ux*ux + uy*uy;
    if ( ll < MMG5_EPSD ) continue;
    li = 1.0 / sqrt(ll);

    /* Tangent vector at p1 */
    if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) {
      t1[0] = li*ux;
      t1[1] = li*uy;
    }
    else {
      t1[0] = p1->n[1];
      t1[1] = -p1->n[0];
    }

    /* Tangent vector at p2 */
    if ( MG_SIN(p2->tag) || p2->tag & MG_NOM ) {
      t2[0] = li*ux;
      t2[1] = li*uy;
    }
    else {
      t2[0] = p2->n[1];
      t2[1] = -p2->n[0];
    }

    /* Calculation of the two Bezier coefficients along the boundary curve */
    ps1   = ux*t1[0] + uy*t1[1];
    b1[0] = p1->c[0] + MMG5_ATHIRD*ps1*t1[0];
    b1[1] = p1->c[1] + MMG5_ATHIRD*ps1*t1[1];

    ps2   = ux*t2[0] + uy*t2[1];
    b2[0] = p2->c[0] - MMG5_ATHIRD*ps2*t2[0];
    b2[1] = p2->c[1] - MMG5_ATHIRD*ps2*t2[1];

    ps1 *= ps1;
    ps2 *= ps2;

    if ( ps1 < MMG5_EPSD || ps2 < MMG5_EPSD ) continue;

    /* \gamma^{\prime\prime}(0); \gamma^\prime(0) = ps*t1 by construction */
    gpp1[0] = 6.0*(p1->c[0] - 2.0*b1[0] + b2[0]);
    gpp1[1] = 6.0*(p1->c[1] - 2.0*b1[1] + b2[1]);

    /* Vector product gpp1 ^ t1 */
    pv = gpp1[0]*t1[1] - gpp1[1]*t1[0];
    M1 = fabs(pv)/ps1;

    /* \gamma^{\prime\prime}(1); \gamma^\prime(1) = -ps*t2 by construction */
    gpp2[0] = 6.0*(p2->c[0] - 2.0*b2[0] + b1[0]);
    gpp2[1] = 6.0*(p2->c[1] - 2.0*b2[1] + b1[1]);

    /* Vector product gpp2 ^ t2 */
    pv = gpp2[0]*t2[1] - gpp2[1]*t2[0];
    M2 = fabs(pv)/ps2;

    M1 = MG_MAX(M1,M2);
    if ( M1 < MMG5_EPSD) continue;
    else {
      ltmp = 8.0*hausd / M1;
      lm = MG_MAX(sqhmin,MG_MIN(ltmp,lm));
    }
  }

  /* Expression of the metric tensor diag(lm,hmax) (in the (tau, n) basis) in the canonical basis */
  n = &p0->n[0];
  sqhmax = 1.0 / sqhmax;
  lm = 1.0 / lm;

  m[0] = lm*n[1]*n[1] + sqhmax*n[0]*n[0];
  m[1] = n[0]*n[1]*(sqhmax-lm);
  m[2] = lm*n[0]*n[0] + sqhmax*n[1]*n[1];

  p0->flag = 2;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 0 if fail, 1 if success
 *
 * Definition of an anisotropic metric tensor field based on the geometry of the
 * domain; this tensor field is intersected by a user-defined tensor field
 */
int MMG2D_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_pPar      ppa;
  double         mm[3],mr[3],isqhmax;
  int            k,l,ip;
  int8_t         ismet;
  char           isdef,i;

  if ( !MMG5_defsiz_startingMessage (mesh,met,__func__) ) {
    return 0;
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
    ppt->s    = 0;
  }

  /* Allocate the structure */
  if ( met->m )
    ismet = 1;
  else {
    ismet = 0;
    if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,3) ) {
      return 0;
    }
  }

  /** Step 1: Set metric at points belonging to a required edge: compute the
   * metric as the mean of the length of the required eges passing through the
   * point */
  if ( !MMG2D_set_metricAtPointsOnReqEdges ( mesh,met,ismet ) ) {
    return 0;
  }

  /* Step 2: Travel all the points (via triangles) in the mesh and set metric tensor */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 ) continue;

    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ppt = &mesh->point[ip];
      if ( !MG_VOK(ppt) || ppt->flag ) continue;
      if ( ismet )
        memcpy(mm,&met->m[3*ip],3*sizeof(double));

      isdef = 0;
      /* Calculation of a metric tensor depending on the anisotropic features of the mesh */
      /* At a singular point, an isotropic metric with size hmax is defined */
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) {
        /* Set the point flag to 1 */
        if ( MMG2D_defaultmet_2d(mesh,met,k,i) ) isdef = 1;
      }
      else if ( MG_EDG(ppt->tag) ) {
        /* Set the point flag to 2 */
        if ( MMG2D_defmetbdy_2d(mesh,met,k,i) ) isdef = 1;
      }

      /* If ppt is an interior point, or if it is a boundary point and the special definition of
       a metric tensor has failed, define a default isotropic metric at ppt */
      if ( !isdef ) {
        MMG2D_defaultmet_2d(mesh,met,k,i);
      }

      /* If a metric is supplied by the user, intersect it with the geometric one */
      if ( ismet && MMG5_intersecmet22(mesh,&met->m[3*ip],mm,mr) )
        memcpy(&met->m[3*ip],mr,3*sizeof(double));
    }
  }

  /** For points with flag 1 (metrec computed by defaultmet_2d), truncation by
   * the local parameters */
  if ( mesh->info.npar ) {
    /* Minimum size feature imposed by triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      /* Retrieve local parameters associated to triangle k */
      for (l=0; l<mesh->info.npar; l++) {
        ppa = &mesh->info.par[l];
        if ( ppa->elt == MMG5_Triangle && ppa->ref == pt->ref ) {
          for (i=0; i<3; i++) {
            ip = pt->v[i];
            if ( mesh->point[ip].flag > 1 ) continue;

            isqhmax = 1./(ppa->hmax*ppa->hmax);
            mm[0] = mm[2] = isqhmax;

            if ( MMG5_intersecmet22(mesh,&met->m[3*ip],mm,mr) ) {
              memcpy(&met->m[3*ip],mr,3*sizeof(double));
            }
          }
          break;
        }
      }
    }
    /* Minimum size feature imposed by vertices */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( (!MG_VOK(ppt)) || ppt->flag > 1 ) continue;

      /* Retrieve local parameters associated to vertex k */
      for (l=0; l<mesh->info.npar; l++) {
        ppa = &mesh->info.par[l];
        if ( ppa->elt == MMG5_Vertex && ppa->ref == ppt->ref ) {
          isqhmax = 1./(ppa->hmax*ppa->hmax);
          mm[0] = mm[2] = isqhmax;

          if ( MMG5_intersecmet22(mesh,&met->m[3*k],mm,mr) ) {
            memcpy(&met->m[3*k],mr,3*sizeof(double));
          }
          break;
        }
      }
    }
  }

  return 1;
}

/**
 * \param dm eigenvalues of the first matrix
 * \param dn eigenvalues of the second matrix
 * \param difsiz maximal size gap authorized by the gradation.
 * \param dir direction in which the sizes are graded.
 * \param ier flag of the modified eigenvalue: (ier & 1) if dm is altered, and (ier & 2) if dn is altered.
 *
 *  Gradation of sizes = 1/sqrt(eigenv of the tensors) in the \a idir direction.
 *
 */
static inline
void MMG2D_gradEigenv(double dm[2],double dn[2],double difsiz,int8_t dir,int8_t *ier) {
  double hm,hn;

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the first direction */
  hm = 1.0 / sqrt(dm[dir]);
  hn = 1.0 / sqrt(dn[dir]);

  if ( hn > hm + difsiz + MMG5_EPSOK ) {
    hn = hm+difsiz;
    dn[dir] = 1.0 / (hn*hn);
    (*ier) = (*ier) | 2;
  }
  else if ( hm > hn + difsiz + MMG5_EPSOK ) {
    hm = hn+difsiz;
    dm[dir] = 1.0 / (hm*hm);
    (*ier) = (*ier) | 1;
  }
}

/**
 * \param dm eigenvalues of the first matrix
 * \param dn eigenvalues of the second matrix
 * \param difsiz maximal size gap authorized by the gradation.
 * \param dir direction in which the sizes are graded.
 * \param ier 2 if dn has been updated, 0 otherwise.
 *
 *  Gradation of size dn = 1/sqrt(eigenv of the tensor) for required points in
 *  the \a idir direction.
 *
 */
static inline
void MMG2D_gradEigenvreq(double dm[2],double dn[2],double difsiz,int8_t dir,int8_t *ier) {
  double hm,hn;

  hm = 1.0 / sqrt(dm[dir]);
  hn = 1.0 / sqrt(dn[dir]);

  if ( hn > hm + difsiz + MMG5_EPSOK ) {
    /* Decrease the size in \a ipslave */
    hn = hm+difsiz;
    dn[dir] = 1.0 / (hn*hn);
    (*ier) = 2;
  }
  else if ( hn + MMG5_EPSOK < hm - difsiz ) {
    /* Increase the size in \a ipslave */
    hn = hm-difsiz;
    dn[dir] = 1.0 / (hn*hn);
    (*ier) = 2;
  }
}

/**
 * \param m first matrix
 * \param n second matrix
 * \param dm eigenvalues of m in the coreduction basis
 * \param dn eigenvalues of n in the coreduction basis
 * \param vp coreduction basis
 * \param ier flag of the updated sizes: (ier & 1) if we dm has been modified, (ier & 2) if dn has been modified.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored in
 * columns
 *
 */
static inline
int MMG2D_updatemet_ani(double *m,double *n,double dm[2],double dn[2],
                         double vp[2][2],int8_t ier ) {
  double det,ip[4];

  det = vp[0][0]*vp[1][1] - vp[0][1]*vp[1][0];
  if ( fabs(det) < MMG5_EPS )  return 0;
  det = 1.0 / det;

  ip[0] =  vp[1][1]*det;
  ip[1] = -vp[1][0]*det;
  ip[2] = -vp[0][1]*det;
  ip[3] =  vp[0][0]*det;

  if ( ier & 1 ) {
    m[0] = dm[0]*ip[0]*ip[0] + dm[1]*ip[2]*ip[2];
    m[1] = dm[0]*ip[0]*ip[1] + dm[1]*ip[2]*ip[3];
    m[2] = dm[0]*ip[1]*ip[1] + dm[1]*ip[3]*ip[3];
  }
  if ( ier & 2 ) {
    n[0] = dn[0]*ip[0]*ip[0] + dn[1]*ip[2]*ip[2];
    n[1] = dn[0]*ip[0]*ip[1] + dn[1]*ip[2]*ip[3];
    n[2] = dn[0]*ip[1]*ip[1] + dn[1]*ip[3]*ip[3];
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param ip1 first edge extremity
 * \param ip2 second edge extremity
 * \param difsiz maximal size gap authorized by the gradation.
 *
 * \return 0 if fail or we don't need to modify the sizes. ier, where (ier & 1)
 * if metric of \a ip1 is altered, and (ier & 2) if metric of \a ip2 is altered.
 *
 * Perform simultaneous reduction of metrics at ip1 points and ip2, and truncate
 * characteristic sizes so that the difference between two corresponding sizes
 * is less than difsiz.
 *
 * Ref : https://www.rocq.inria.fr/gamma/Frederic.Alauzet/cours/cea2010_V2.pdf
 *
 */
int MMG2D_grad2met_ani(MMG5_pMesh mesh,MMG5_pSol met,int ip1,int ip2,double difsiz) {
  double       dm[2],dn[2];
  double       vp[2][2],*m,*n;
  int8_t       ier;

  ier = 0;

  m = &met->m[met->size*ip1];
  n = &met->m[met->size*ip2];

  /* Simultaneous reduction of m1 and m2 */
  if ( !MMG2D_simred(mesh,m,n,dm,dn,vp) ) {
    return 0;
  }

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the first direction */
  MMG2D_gradEigenv(dm,dn,difsiz,0,&ier);

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the second direction */
  MMG2D_gradEigenv(dm,dn,difsiz,1,&ier);

  if ( !ier ) {
    return 0;
  }

  /* Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored in
   * columns */
  if ( !MMG2D_updatemet_ani(m,n,dm,dn,vp,ier ) ) {
    return 0;
  }

  return ier;

}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param ipmaster edge extremity that cannot be modified
 * \param ipslave edge extremity to modify to respect the gradation.
 * \param difsiz maximal size gap authorized by the gradation.
 *
 * \return 0 if fail or we don't need to update the size of \a ipslave, 1 if
 * its size has been updated.
 *
 * Perform simultaneous reduction of metrics at ipmaster points and ipslave, and
 * modify the characteristic size of \a ipslave so that the difference between
 * the two sizes is less than difsiz.
 *
 * Ref : https://www.rocq.inria.fr/gamma/Frederic.Alauzet/cours/cea2010_V2.pdf
 *
 */
static inline
int MMG2D_grad2metreq_ani(MMG5_pMesh mesh,MMG5_pSol met,int ipmaster,int ipslave,double difsiz) {
  double       dm[2],dn[2];
  double       vp[2][2],*m,*n;
  int8_t       ier;

  ier = 0;

  m = &met->m[met->size*ipmaster];
  n = &met->m[met->size*ipslave];

  /* Simultaneous reduction of m1 and m2 */
  if ( !MMG2D_simred(mesh,m,n,dm,dn,vp) ) {
    return 0;
  }

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the first direction */
  MMG2D_gradEigenvreq(dm,dn,difsiz,0,&ier);

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the second direction */
  MMG2D_gradEigenvreq(dm,dn,difsiz,1,&ier);

  if ( !ier ) {
    return 0;
  }

  /* Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored in
   * columns */
  if ( !MMG2D_updatemet_ani(m,n,dm,dn,vp,ier ) ) {
    return 0;
  }

  return ier;

}


/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * Anisotropic mesh gradation routine
 *
 */
int MMG2D_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            hgrad,ll,difsiz;
  int               k,it,ip1,ip2,maxit,nup,nu;
  char              i,i1,i2,ier;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n");

  /** Mark the edges belonging to a required entity */
  MMG5_mark_pointsOnReqEdge_fromTria ( mesh );

  for (k=1; k<=mesh->np; k++) {
    mesh->point[k].flag = mesh->base;
  }

  hgrad = mesh->info.hgrad;
  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];

        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;

        /* Skip points belonging to a required edge */
        if ( p1->s || p2->s ) continue;

        ll = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0])
          + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);
        ll = sqrt(ll);

        /* Maximum allowed difference between the prescribed sizes in p1 and p2 */
        difsiz = ll*hgrad;

        /* bit 0 of ier = 0 if metric at point ip1 is untouched, 1 otherwise;
         * bit 1 of ier = 0 if metric at point ip2 is untouched, 1 otherwise */
        ier = MMG2D_grad2met_ani(mesh,met,ip1,ip2,difsiz);
        if ( ier & 1 ) {
          p1->flag = mesh->base;
          nu++;
        }
        if ( ier & 2 ) {
          p2->flag = mesh->base;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  }

  return 1;

}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return the number of updated metrics
 *
 * Anisotropic gradation on required points.
 *
 */
int MMG2D_gradsizreq_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            hgrad,ll,difsiz;
  int               k,it,ip1,ip2,ipslave,ipmaster,maxit,nup,nu;
  int8_t            ier;
  char              i,i1,i2;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading required points.\n");
  }

  if ( mesh->info.hgrad < 0. ) {
    /** Mark the edges belonging to a required entity */
    MMG5_mark_pointsOnReqEdge_fromTria ( mesh );
  }

  hgrad = mesh->info.hgrad;
  it = nup = 0;
  maxit = 100;

  do {
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];

        if ( abs ( p1->s - p2->s ) < 2 ) {
          /* No size to propagate */
          continue;
        }
        else if ( p1->s > p2->s ) {
          ipmaster = ip1;
          ipslave  = ip2;
        }
        else {
          assert ( p2->s > p1->s );
          ipmaster = ip2;
          ipslave  = ip1;
        }

        ll = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0])
          + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);
        ll = sqrt(ll);

        /* Maximum allowed difference between the prescribed sizes in p1 and p2 */
        difsiz = ll*hgrad;

        /* Impose the gradation to ipslave from ipmaster */
        ier = MMG2D_grad2metreq_ani(mesh,met,ipmaster,ipslave,difsiz);

        if ( ier ) {
          mesh->point[ipslave].s = mesh->point[ipmaster].s - 1;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 && nup ) {
    fprintf(stdout,"     gradation (required): %7d updated, %d iter.\n",nup,it);
  }
  return nup;
}
