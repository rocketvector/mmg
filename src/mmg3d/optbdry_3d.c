/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * \file mmg3d/optbdry_3d.c
 * \brief Functions for the optimization of very bad elements.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param octree pointer toward the octree
 * \param k index of the tetra of bad quality that we try to improve
 *
 * \return -1 if fail, 1 if we move a point, 0 if not.
 *
 * Try to improve a tetra of very bad quality by moving its vertices.
 *
 */
int MMG3D_movetetrapoints(MMG5_pMesh mesh,MMG5_pSol met,_MMG3D_pOctree octree,int k) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  double        *n;
  int           i,j,i0,ier,lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv;
  int           improve,internal,nm,maxit,base,ns;

  improve = 1;
  internal = 1;
  nm = ns = 0;
  maxit = 1;
  base = mesh->base;

  pt = &mesh->tetra[k];

  /* point j on face i */
  for (i=0; i<4; i++) {
    for (j=0; j<3; j++) {
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pxt->tag[_MMG5_iarf[i][j]] & MG_REQ )  continue;
      }
      else  pxt = 0;
      i0  = _MMG5_idir[i][j];
      ppt = &mesh->point[pt->v[i0]];
      if ( ppt->flag == base )  continue;
      else if ( MG_SIN(ppt->tag) )  continue;

      if ( maxit != 1 ) {
        ppt->flag = base;
      }
      ier = 0;
      if ( ppt->tag & MG_BDY ) {
        /* Catch a boundary point by a boundary face */
        if ( !pt->xt || !(MG_BDY & pxt->ftag[i]) )  continue;
        else if( ppt->tag & MG_NOM ){
          if( mesh->adja[4*(k-1)+1+i] ) continue;
          ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,1);
          if( !ier )  continue;
          else if ( ier>0 )
            ier = _MMG5_movbdynompt(mesh,met,octree,listv,ilistv,lists,ilists,improve);
          else
            return(-1);
        }
        else if ( ppt->tag & MG_GEO ) {
          ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
          if ( !ier )  continue;
          else if ( ier>0 )
            ier = _MMG5_movbdyridpt(mesh,met,octree,listv,ilistv,lists,ilists,improve);
          else
            return(-1);
        }
        else if ( ppt->tag & MG_REF ) {
          ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
          if ( !ier )
            continue;
          else if ( ier>0 )
            ier = _MMG5_movbdyrefpt(mesh,met,octree,listv,ilistv,lists,ilists,improve);
          else
            return(-1);
        }
        else {
          ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
          if ( !ier )
            continue;
          else if ( ier<0 )
            return(-1);

          n = &(mesh->xpoint[ppt->xp].n1[0]);

          if ( !MG_GET(pxt->ori,i) ) {
            if ( !_MMG5_directsurfball(mesh,pt->v[i0],lists,ilists,n) )
              continue;
          }
          ier = _MMG5_movbdyregpt(mesh,met, octree, listv,ilistv,lists,ilists,improve);
          if ( ier )  ns++;
        }
      }
      else if ( internal ) {
        ilistv = _MMG5_boulevolp(mesh,k,i0,listv);
        if ( !ilistv )  continue;
        ier = _MMG5_movintpt(mesh,met,octree,listv,ilistv,improve);
      }
      if ( ier ) {
        nm++;
        if(maxit==1){
          ppt->flag = base;
        }
      }
    }
  }

  return ( nm > 0 );
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param k index of the tetra of bad quality that we try to improve
 * \param i index of point to collapse
 *
 * \return -1 if fail, 1 if we collapse a point, 0 if not.
 *
 * Try to remove point i of tet k, try the three edges of k containing i
 *
 */
int _MMG3D_coledges(MMG5_pMesh mesh,MMG5_pSol met,int k,int i) {
  MMG5_pTetra pt;
  double      len;
  int         ied,iedg,iq,i1,ilistcol,listcol[MMG3D_LMAX+2];
  int         ier;
  char        iface,ief;

  pt = &mesh->tetra[k];

   /*3 possibilities to remove the vertex ib*/
    for(ied = 0 ; ied<3 ;ied++) {
      iedg  = _MMG5_arpt[i][ied];
      len =  _MMG5_lenedg(mesh,met,iedg,pt);
#warning unable to check if the edge is req without computing the shell because the tet may be non-boundary but have a boundary edge
      if(len > 1.1) continue;
      iface = _MMG5_ifar[iedg][0];
      ief   = _MMG5_iarfinv[iface][iedg];
      iq    = _MMG5_idir[iface][_MMG5_iprv2[ief]];
      if(iq==i) {
        iface = _MMG5_ifar[iedg][1];
        ief   = _MMG5_iarfinv[iface][iedg];
        iq    = _MMG5_idir[iface][_MMG5_iprv2[ief]];
      }
      i1    = _MMG5_idir[iface][_MMG5_inxt2[ief]];

      ilistcol = _MMG5_boulevolp(mesh,k,i1,listcol);

      ilistcol = _MMG5_chkcol_int(mesh,met,k,iface,ief,listcol,ilistcol,2);
      if ( ilistcol > 0 ) {
        ier = _MMG5_colver(mesh,met,listcol,ilistcol,iq,2);
        if ( ilistcol < 0 ) continue;
        if ( ier < 0 ) return(-1);
        else if(ier) {
          _MMG3D_delPt(mesh,ier);
          //printf("del succeed\n");
          //exit(0);
          return(1);
        }
      }
    }
    return(0);
}


/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param octree pointer toward the octree
 * \param k index of the tetra of bad quality that we try to improve
 * \param i index of the point that we want to delete.
 *
 * \return -1 if fail, 1 if we move a point, 0 if not.
 *
 * Try to delete the point \a i of tetra \a k. Try all the edges containing \a i.
 *
 */
int _MMG3D_deletePoint(MMG5_pMesh mesh,  MMG5_pSol met,_MMG3D_pOctree octree,
                       int k,int i) {
  int         il,ilist,iel,ip,list[MMG3D_LMAX+2];

  ilist = _MMG5_boulevolp(mesh,k,i,list);

  for(il = 0 ; il<ilist ; il++) {
    iel = list[il] / 4;
    ip  = list[il] % 4;
    if( _MMG3D_coledges(mesh,met,iel,ip) > 0 ) {
      return(1);
    }
  }

  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param k   index of a tetra
 *
 * \return -1 if fail, 0 if we do not improve the tetra, 1 otherwise.
 *
 * Try to optimize the tetra k. This tetra has a face on the boundary.
 *
 */
int MMG3D_optbdry(MMG5_pMesh mesh,MMG5_pSol met,_MMG3D_pOctree octree,int k) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  double       lint,lbdy,len,lmin,lmax;
  int          imax;
  int          j,list[MMG3D_LMAX+2];
  int          iedg,ier,ilist,ied,ia,it1,it2,ret,improved;
  char         ifac2,ifac;

  pt = &mesh->tetra[k];
  assert(pt->xt);

  pxt = &mesh->xtetra[pt->xt];

  for(ifac=0 ; ifac<4 ; ifac++)
    if ( pxt->ftag[ifac] & MG_BDY ) break;

  assert ( ifac!=4 );

  ier = 0;

  /*First : try to move the vertex in order to improve the quality*/
  improved = 0;
  if ( !mesh->info.noinsert ) {
    for(j = 0 ; j<10 ; j++) {
      ier = MMG3D_movetetrapoints(mesh,met,octree,k);
      if ( ier<0 ) return -1;
      else if ( !ier ) break;
      improved += ier;
    }

    if ( improved ) {
      if ( mesh->info.ddebug || mesh->info.imprim > 9 ) {
        printf(" -- optbdry: %d move in tetra %d\n",improved,k);
      }
    }
  }

  if ( !mesh->info.noinsert ) {
    /*look at the lenght edges*/
    imax = -1;
    lmax = 0;
    lint = 0;
#warning how to properly detect the boundary edges
    for ( ifac2 = 0 ; ifac2<4 ;ifac2++) {
      if ( pxt->ftag[ifac2] & MG_BDY ) continue;

      for ( ied = 0; ied < 3; ++ied ) {
#warning this can not work because a tetra may have multiple boundary faces
        iedg  = _MMG5_arpt[ifac][ied];
      len =  _MMG5_lenedg(mesh,met,iedg,pt);
      lint += len;
      if(len > lmax) {
        imax = iedg;
        lmax = len;
      }
    }
    lint /= 3.;
    lbdy = 0;
    lmin = 1000;
    for(ied = 0 ; ied<3 ;ied++) {
      iedg  = _MMG5_iarf[ifac][ied];
      len = _MMG5_lenedg(mesh,met,iedg,pt);
      lbdy += len;
      if(len < lmin) {
        lmin = len;
      }
    }
    lbdy /= 3.;

    /** Try to split too long internal edges */
    if(lint > 2.*lbdy) {
      ier = _MMG3D_splitItem(mesh,met,octree,k,imax,1.00001);
      if ( ier ) return 1;

      for(ied = 0 ; ied<3 ;ied++) {
        iedg  = _MMG5_arpt[ifac][ied];
        if(iedg == imax) continue;
        ier = _MMG3D_splitItem(mesh,met,octree,k,iedg,1.00001);
        if(ier)  {
          return 1;
        }
      }

      /* try to remove the non-bdry vertex : with all the edges containing the
       * vertex*/
#warning can not work because the edge may be boundary
      ier = 0;//_MMG3D_deletePoint(mesh,met,octree,k,i);
      if ( ier < 0 )  return -1;
      else if ( ier ) return  1;
    }
  }

  /*First : try to swap the 3 other edges*/
  if(!mesh->info.noswap) {
    for(ied = 0 ; ied<3 ;ied++) {
      iedg  = _MMG5_arpt[ifac][ied];
      ier = _MMG3D_swpItem(mesh,met,octree,k,iedg);
      if ( ier < -1 ) return -1;
      else if ( ier ) return 1;
      }
    }

    /*try to swap the 3 face edges*/
    for (j=0; j<3; j++) {
      ia  = _MMG5_iarf[ifac][j];

      /* No swap of geometric edge */
      if ( MG_EDG(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) ||
           (pxt->tag[ia] & MG_NOM) )
        continue;

      ret = _MMG5_coquilface(mesh,k,ifac,ia,list,&it1,&it2,0);
      ilist = ret / 2;
      if ( ret < 0 )  return(-1);
      /* CAUTION: trigger collapse with 2 elements */
      if ( ilist <= 1 )  continue;
      ier = _MMG5_chkswpbdy(mesh,met,list,ilist,it1,it2,2);
      if ( ier <  0 )
        return -1;
      else if ( ier ) {
        ier = _MMG5_swpbdy(mesh,met,list,ret,it1,octree,2);
        if ( ier < 0 )  return(-1);
        else if(ier) {
          return(1);
        }
      }
    }
  }

  /*Try : try to remove the non-bdry vertex*/
  if ( !mesh->info.noinsert ) {
    ier = 0;//_MMG3D_coledges(mesh,met,k,ifac);

    if ( !ier ) {
      /* try to remove the non-bdry vertex : with all the edges containing the
       * vertex*/
      ier = 0;//_MMG3D_deletePoint(mesh,met,octree,k,ifac);
    }
  }

  return ( ier || improved );
}
