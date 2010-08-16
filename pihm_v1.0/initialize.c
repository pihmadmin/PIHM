/*******************************************************************************
 * File        : initialize.c                                                  *
 * Function    : initialization of data structue                               *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               * 
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * Part of initialization work is done in read_alloc, and the rest is done here*
 *                                                                             *
 * This code is free for users with research purpose only, if appropriate      *
 * citation is refered. However, there is no warranty in any format for this   *
 * product.                                                                    *
 *                                                                             *
 * For questions or comments, please contact the authors of the reference.     *
 * One who want to use it for other consideration may also contact Dr.Duffy    *
 * at cxd11@psu.edu.                                                           *    
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sundialstypes.h"
#include "nvector_serial.h"
#include "ihm.h"  

void initialize(char *filename, Model_Data DS, Control_Data *CS, N_Vector CV_Y)
{
  int i;
  realtype a_x, a_y, b_x, b_y, c_x, c_y;
  realtype a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax; 
  realtype tempvalue;
  
  FILE *int_file;
  char *fn;
  
  printf("\nInitializing data structure ... ");
  
  for(i=0; i<DS->NumEle; i++)
  {
    a_x = DS->Node[DS->Ele[i].node[0]-1].x;
    b_x = DS->Node[DS->Ele[i].node[1]-1].x;
    c_x = DS->Node[DS->Ele[i].node[2]-1].x;
    a_y = DS->Node[DS->Ele[i].node[0]-1].y;
    b_y = DS->Node[DS->Ele[i].node[1]-1].y;
    c_y = DS->Node[DS->Ele[i].node[2]-1].y;
    a_zmin = DS->Node[DS->Ele[i].node[0]-1].zmin;
    b_zmin = DS->Node[DS->Ele[i].node[1]-1].zmin;
    c_zmin = DS->Node[DS->Ele[i].node[2]-1].zmin;
    a_zmax = DS->Node[DS->Ele[i].node[0]-1].zmax;
    b_zmax = DS->Node[DS->Ele[i].node[1]-1].zmax;
    c_zmax = DS->Node[DS->Ele[i].node[2]-1].zmax;
    
    DS->Ele[i].area = 0.5*((b_x - a_x)*(c_y - a_y) - (b_y - a_y)*(c_x - a_x));
    DS->Ele[i].zmin = (a_zmin + b_zmin + c_zmin)/3.0;
    DS->Ele[i].zmax = (a_zmax + b_zmax + c_zmax)/3.0;
    
    DS->Ele[i].edge[0] = pow((a_x - b_x), 2) + pow((a_y - b_y), 2);
    DS->Ele[i].edge[1] = pow((b_x - c_x), 2) + pow((b_y - c_y), 2);
    DS->Ele[i].edge[2] = pow((c_x - a_x), 2) + pow((c_y - a_y), 2);
    
    /* calculate centroid of triangle */
    /*     
    DS->Ele[i].x = (a_x + b_x + c_x)/3.0;
    DS->Ele[i].y = (a_y + b_y + c_y)/3.0;
    */
    
    /* calculate circumcenter of triangle */
    DS->Ele[i].x = a_x - ((b_y - a_y)*DS->Ele[i].edge[2] - (c_y - a_y)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
    DS->Ele[i].y = a_y + ((b_x - a_x)*DS->Ele[i].edge[2] - (c_x - a_x)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
    
    DS->Ele[i].edge[0] = sqrt(DS->Ele[i].edge[0]);
    DS->Ele[i].edge[1] = sqrt(DS->Ele[i].edge[1]);
    DS->Ele[i].edge[2] = sqrt(DS->Ele[i].edge[2]);
  }
  
  /* allocate flux */
  DS->FluxSurf = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  DS->FluxSub = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  DS->FluxRiv = (realtype **)malloc(DS->NumRiv*sizeof(realtype));
  DS->EleET = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  
  for(i=0; i<DS->NumEle; i++)
  {
    DS->FluxSurf[i] = (realtype *)malloc(3*sizeof(realtype));
    DS->FluxSub[i] = (realtype *)malloc(3*sizeof(realtype));
    DS->EleET[i] = (realtype *)malloc(4*sizeof(realtype));
  }
  
  for(i=0; i<DS->NumRiv; i++)
  {
    DS->FluxRiv[i] = (realtype *)malloc(6*sizeof(realtype));
  }
  
  DS->ElePrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->EleVic = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->Recharge = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->EleIS = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->EleISmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->EleETP = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->Ele2IS = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  DS->EleNetPrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
     
  for(i=0; i<DS->NumEle; i++)
  {
    DS->Ele[i].Ksat = DS->Soil[(DS->Ele[i].soil-1)].Ksat;
    DS->Ele[i].Porosity = DS->Soil[(DS->Ele[i].soil-1)].SitaS - 
                          DS->Soil[(DS->Ele[i].soil-1)].SitaR;
    DS->Ele[i].Alpha = DS->Soil[(DS->Ele[i].soil-1)].Alpha;
    DS->Ele[i].Beta = DS->Soil[(DS->Ele[i].soil-1)].Beta; 
    DS->Ele[i].Sf = DS->Soil[(DS->Ele[i].soil-1)].Sf;
    DS->Ele[i].Rough = DS->Soil[DS->Ele[i].soil-1].Rough;
  }
  
  for(i=0; i<DS->NumRiv; i++)
  {
    DS->Riv[i].x = (DS->Node[DS->Riv[i].FromNode-1].x + 
                    DS->Node[DS->Riv[i].ToNode-1].x)/2;
    DS->Riv[i].y = (DS->Node[DS->Riv[i].FromNode-1].y + 
                    DS->Node[DS->Riv[i].ToNode-1].y)/2;
    DS->Riv[i].zmax = (DS->Node[DS->Riv[i].FromNode-1].zmax + DS->Node[DS->Riv[i].ToNode-1].zmax)/2;
    DS->Riv[i].depth = DS->Riv_Shape[DS->Riv[i].shape-1].depth;
    DS->Riv[i].zmin = DS->Riv[i].zmax - DS->Riv[i].depth;
    
    DS->Riv[i].Length = sqrt(pow(DS->Node[DS->Riv[i].FromNode-1].x - 
                                 DS->Node[DS->Riv[i].ToNode-1].x, 2) + 
                             pow(DS->Node[DS->Riv[i].FromNode-1].y - 
                                 DS->Node[DS->Riv[i].ToNode-1].y, 2));
  }
  
  /* initialize state varible */
  /* relax cases */
  if (CS->int_type == 0)
  {
    for(i=0; i<DS->NumEle; i++)
    {
      DS->EleIS[i] = 0;
      NV_Ith_S(CV_Y, i) = 0;
      NV_Ith_S(CV_Y, i + DS->NumEle) = 0.08;
      NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele[i].zmax - DS->Ele[i].zmin - 0.1;
    }  
    
    for(i=0; i<DS->NumRiv; i++)
    {
      NV_Ith_S(CV_Y, i + 3*DS->NumEle) = 0;
    }
  }
  /* type mode */  
  else if (CS->int_type == 1)
  {
    for(i=0; i<DS->NumEle; i++)
    {
      DS->EleIS[i] = DS->Ele_IC[DS->Ele[i].IC-1].interception;
      NV_Ith_S(CV_Y, i) = DS->Ele_IC[DS->Ele[i].IC-1].surf;
      NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[DS->Ele[i].IC-1].unsat;
      NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele_IC[DS->Ele[i].IC-1].sat;
      
      if ((NV_Ith_S(CV_Y, i + DS->NumEle) + NV_Ith_S(CV_Y, i + 2*DS->NumEle)) >= (DS->Ele[i].zmax - DS->Ele[i].zmin))
      {
        NV_Ith_S(CV_Y, i + DS->NumEle) = ((DS->Ele[i].zmax - DS->Ele[i].zmin) - NV_Ith_S(CV_Y, i + 2*DS->NumEle))*0.9;
        if (NV_Ith_S(CV_Y, i + DS->NumEle) < 0) {NV_Ith_S(CV_Y, i + DS->NumEle) = 0; }
      } 
    }  
    
    for(i=0; i<DS->NumRiv; i++)
    {
      NV_Ith_S(CV_Y, i + 3*DS->NumEle) = DS->Riv_IC[DS->Riv[i].IC-1].value;
    }
  }  
  /* hot start mode */
  else
  {
    fn = (char *)malloc((strlen(filename)+4)*sizeof(char));
    strcpy(fn, filename);
    int_file = fopen(strcat(fn, ".int"), "r");
  
    if(int_file == NULL)
    {
      printf("\n  Fatal Error: %s.int is in use or does not exist!\n", filename);
      exit(1);
    }
    else
    {
      for(i=0; i<DS->NumEle; i++)
      {
        fscanf(int_file, "%lf", &tempvalue);
        if(tempvalue <= 0) {tempvalue = 0.01;}
        NV_Ith_S(CV_Y, i + DS->NumEle) = tempvalue;
      }
      
      for(i=0; i<DS->NumEle; i++)
      {
        fscanf(int_file, "%lf", &tempvalue);
        if(tempvalue <= 0) {tempvalue = 0.01;}
        if(tempvalue >= (DS->Ele[i].zmax - DS->Ele[i].zmin)) {tempvalue = (DS->Ele[i].zmax - DS->Ele[i].zmin) - 0.01;}
        NV_Ith_S(CV_Y, i + 2*DS->NumEle) = tempvalue;
      } 
      
      for(i=0; i<DS->NumEle; i++)
      {
        DS->EleIS[i] = 0;
        NV_Ith_S(CV_Y, i) = 0;
      }  
    
      for(i=0; i<DS->NumRiv; i++)
      {
        NV_Ith_S(CV_Y, i + 3*DS->NumEle) = 0;
      } 
    }
    fclose(int_file); 
  }
  printf("done.\n");
}

