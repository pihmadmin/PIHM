/*******************************************************************************
 * File        : print.c	                                               *
 * Version     : Nov, 2007 (2.0)                                               *
 * Function    : print out model results output files                          *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) This file is downgraded from Version 1.0, as no ancillary results are    *
 *    being output			                                       *
 * b) Only state variables and flux to/in/accross river and its bed are being  *
 *    output							               *
 * c) Addition of Average Function to output average variables at regular time *
 *    intervals								       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"  
#include "cvode.h" 
#include "cvdense.h"
/*Temporal average of State vectors */
void avgResults_NV(FILE *fpin,realtype *tmpVarCal,N_Vector tmpNV,int tmpIntv, int tmpNumObj,realtype tmpt,int tmpInitObj)
        {
        int j;
        for(j=0;j<tmpNumObj;j++)
                {
                tmpVarCal[j]=tmpVarCal[j]+NV_Ith_S(tmpNV,j+tmpInitObj);
                }
        if(((int)tmpt%tmpIntv)==0)      
                {
                fprintf(fpin,"%lf\t",tmpt);
                for(j=0;j<tmpNumObj;j++)
                        {               
                        fprintf(fpin,"%lf\t",tmpVarCal[j]/tmpIntv);
                        tmpVarCal[j]=0; 
                        }
                fprintf(fpin,"\n");     
                fflush(fpin);           
                }
        }
/* Temporal average of Derived states */
void avgResults_MD(FILE *fpin,realtype *tmpVarCal,Model_Data tmpDS,int tmpIntv, int tmpNumObj,realtype tmpt,int tmpFC)
        {
        int j;
        switch(tmpFC)
                {
                case 3: 
                case 4:
                case 5:                 
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleET[j][tmpFC-3];
                                }
                        break;
                case 6: 
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleIS[j];
                                }
                        break;
                case 7:
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+(tmpDS->EleSnowCanopy[j]+tmpDS->EleSnowGrnd[j]);
                                }
                        break;
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->FluxRiv[j][tmpFC-8];
                                }
                        break;
                case 19:
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->Recharge[j];
                                }
                        break;
                default:
                        break;
                }
        if(((int)tmpt%tmpIntv)==0)
                {
                fprintf(fpin,"%lf\t",tmpt);
                for(j=0;j<tmpNumObj;j++)
                        {
                        fprintf(fpin,"%lf\t",tmpVarCal[j]/tmpIntv);
                        tmpVarCal[j]=0;
                        }
                fprintf(fpin,"\n");
                fflush(fpin);
                }
        }
/* print individual states */
void PrintData(FILE **outp,Control_Data *cD, Model_Data DS, N_Vector CV_Y, realtype t)
	{
	int k;
	if(cD->gwD==1)
		{
        	avgResults_NV(outp[0],DS->PrintVar[0],CV_Y,cD->gwDInt,DS->NumEle+DS->NumRiv,t,2*DS->NumEle);
		}
        if(cD->surfD==1)
                {
        	avgResults_NV(outp[1],DS->PrintVar[1],CV_Y,cD->surfDInt,DS->NumEle,t,0*DS->NumEle);
		}
	for(k=0;k<3;k++)
		{
		if(cD->et[k]==1)
			{
			avgResults_MD(outp[2+k],DS->PrintVar[2+k],DS,cD->etInt,DS->NumEle,t,3+k);
			}
		}
        if(cD->IsD==1)
                {
        	avgResults_MD(outp[5],DS->PrintVar[5],DS,cD->IsDInt,DS->NumEle,t,6);
                }
        if(cD->snowD==1)
                {
        	avgResults_MD(outp[6],DS->PrintVar[6],DS,cD->snowDInt,DS->NumEle,t,7);
                }
        for(k=0;k<=10;k++)
        	{
		if(cD->rivFlx[k]==1)
			{
        		avgResults_MD(outp[7+k],DS->PrintVar[k+7],DS,cD->rivFlxInt,DS->NumRiv,t,k+8);
			}
        	}
        if(cD->rivStg==1)
                {
        	avgResults_NV(outp[18],DS->PrintVar[18],CV_Y,cD->rivStgInt,DS->NumRiv,t,3*DS->NumEle);
                }
        if(cD->Rech==1)
                {
       		avgResults_MD(outp[20],DS->PrintVar[20],DS,cD->RechInt,DS->NumEle,t,19);
                }
        if(cD->usD==1)
                {
        	avgResults_NV(outp[19],DS->PrintVar[19],CV_Y,cD->usDInt,DS->NumEle,t,1*DS->NumEle);
                }
	
	}
  
