/*********************************************************************************
 * File        : initialize.c                                                    *
 * Function    : initialization of elemental attributes using relational database*
 * Version     : Nov, 2007 (2.0)                                                 *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                    *
 * Developer of PIHM1.0:        Yizhong Qu                                       *
 *-------------------------------------------------------------------------------*
 *                                                                               *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0..............................*
 * a) Correction of edge calculation term                                        *
 * b) Initialization of new variables for new process, shape representations  and*
 *    calibration (e.g. ET, Infiltration, Macropore, Stormflow, Element beneath  *
 *    river, river shapes, river bed property, thresholds for root zone,         *
 *    infiltration and macropore depths, land cover attributes etc)              *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                         *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sundialstypes.h"
#include "nvector_serial.h"
#include "pihm.h"  


void initialize(char *filename, Model_Data DS, Control_Data *CS, N_Vector CV_Y)
	{
  	int i,j,tmpBool,BoolBR,BoolR=0;
  	realtype a_x, a_y, b_x, b_y, c_x, c_y,distX,distY;
  	realtype a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax; 
  	realtype tempvalue1,tempvalue2,tempvalue3;
  	FILE *init_file;
  	char *fn;
  	realtype *zmin_cor;

  	zmin_cor=(realtype *)malloc(DS->NumEle*sizeof(realtype));
  
  	printf("\nInitializing data structure ... ");

  	/* allocate memory storage to flux terms */
  	DS->FluxSurf = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->FluxSub = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->FluxRiv = (realtype **)malloc(DS->NumRiv*sizeof(realtype));
  	DS->EleET = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->ElePrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleViR = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->Recharge = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleIS = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleISmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleISsnowmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleSnow = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleSnowGrnd = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleSnowCanopy = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleTF = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleETloss = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleNetPrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
          
           
	 for(i=0;i<DS->NumEle;i++)
		{
    		DS->FluxSurf[i] = (realtype *)malloc(3*sizeof(realtype));
    		DS->FluxSub[i] = (realtype *)malloc(3*sizeof(realtype));
    		DS->EleET[i] = (realtype *)malloc(4*sizeof(realtype));

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
    		DS->Ele[i].zmax = (a_zmax + b_zmax + c_zmax)/3.0;
    		DS->Ele[i].zmin = (a_zmin + b_zmin + c_zmin)/3.0; 
    		DS->Ele[i].edge[0] = pow((b_x - c_x), 2) + pow((b_y - c_y), 2);
    		DS->Ele[i].edge[1] = pow((c_x - a_x), 2) + pow((c_y - a_y), 2);
    		DS->Ele[i].edge[2] = pow((a_x - b_x), 2) + pow((a_y - b_y), 2);

    		/* calculate centroid of triangle */
    		DS->Ele[i].x = (a_x + b_x + c_x)/3.0;
    		DS->Ele[i].y = (a_y + b_y + c_y)/3.0;
    
    
    		/* calculate circumcenter of triangle */
  		/*  DS->Ele[i].x = a_x - ((b_y - a_y)*DS->Ele[i].edge[2] - (c_y - a_y)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
    		DS->Ele[i].y = a_y + ((b_x - a_x)*DS->Ele[i].edge[2] - (c_x - a_x)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area); 
    		*/
    		DS->Ele[i].edge[0] = sqrt(DS->Ele[i].edge[0]);
    		DS->Ele[i].edge[1] = sqrt(DS->Ele[i].edge[1]);
    		DS->Ele[i].edge[2] = sqrt(DS->Ele[i].edge[2]);
    		DS->Ele[i].KsatH = CS->Cal.KsatH*DS->Geol[(DS->Ele[i].geol-1)].KsatH;
    		DS->Ele[i].KsatV = CS->Cal.KsatV*DS->Geol[(DS->Ele[i].geol-1)].KsatV;
    		DS->Ele[i].infKsatV = CS->Cal.infKsatV*DS->Soil[(DS->Ele[i].soil-1)].KsatV;
    		DS->Ele[i].Porosity = CS->Cal.Porosity*(DS->Soil[(DS->Ele[i].soil-1)].ThetaS - DS->Soil[(DS->Ele[i].soil-1)].ThetaR);
		/* Note above porosity statement should be replaced by geologic porosity (in comments below) if the data is available */
//		DS->Ele[i].Porosity = CS->Cal.Porosity*(DS->Geol[(DS->Ele[i].geol-1)].ThetaS - DS->Geol[(DS->Ele[i].geol-1)].ThetaR);
 		if((DS->Ele[i].Porosity>1)&&(DS->Ele[i].Porosity==0))
			{
			printf("Warning: Porosity value out of bounds");
			getchar();
			} 
		DS->Ele[i].Alpha = CS->Cal.Alpha*DS->Soil[(DS->Ele[i].soil-1)].Alpha;
    		DS->Ele[i].Beta = CS->Cal.Beta*DS->Soil[(DS->Ele[i].soil-1)].Beta; 
		/* Note above van genuchten statement should be replaced by geologic parameters (in comments below) if the data is available */
//		DS->Ele[i].Alpha = CS->Cal.Alpha*DS->Geol[(DS->Ele[i].geol-1)].Alpha;
//    		DS->Ele[i].Beta = CS->Cal.Beta*DS->Geol[(DS->Ele[i].geol-1)].Beta; 
    		DS->Ele[i].hAreaF = CS->Cal.hAreaF*DS->Soil[(DS->Ele[i].soil-1)].hAreaF; 
    		DS->Ele[i].vAreaF = CS->Cal.vAreaF*DS->Geol[(DS->Ele[i].geol-1)].vAreaF; 
    		DS->Ele[i].macKsatV = CS->Cal.macKsatV*DS->Soil[(DS->Ele[i].soil-1)].macKsatV; 
    		DS->Ele[i].macKsatH = CS->Cal.macKsatH*DS->Geol[(DS->Ele[i].geol-1)].macKsatH; 
		DS->Ele[i].macD=CS->Cal.macD*DS->Geol[DS->Ele[i].geol-1].macD;
    		DS->Ele[i].infD=CS->Cal.infD*DS->Soil[DS->Ele[i].soil-1].infD;
    
    		DS->Ele[i].RzD=CS->Cal.RzD*DS->LandC[DS->Ele[i].LC-1].RzD;  
    		DS->Ele[i].LAImax = DS->LandC[DS->Ele[i].LC-1].LAImax;
    		DS->Ele[i].Rmin = DS->LandC[DS->Ele[i].LC-1].Rmin;
    		DS->Ele[i].Rs_ref = DS->LandC[DS->Ele[i].LC-1].Rs_ref;
    		DS->Ele[i].Albedo = CS->Cal.Albedo*DS->LandC[DS->Ele[i].LC-1].Albedo;
		if(DS->Ele[i].Albedo>1)
			{
                        printf("Warning: Albedo out of bounds");
                        getchar();
			}
    		DS->Ele[i].VegFrac = CS->Cal.VegFrac*DS->LandC[DS->Ele[i].LC-1].VegFrac;                
    		DS->Ele[i].Rough = CS->Cal.Rough*DS->LandC[DS->Ele[i].LC-1].Rough;

		DS->Ele[i].windH=DS->windH[DS->Ele[i].WindVel-1];
		}
  	for(i=0; i<DS->NumRiv; i++)
  		{
    		DS->FluxRiv[i] = (realtype *)malloc(11*sizeof(realtype));
		for(j=0;j<3;j++)
			{
			/* Note: Strategy to use BC < -4 for river identification */
			if(DS->Ele[DS->Riv[i].LeftEle-1].nabr[j]==DS->Riv[i].RightEle)
				{
				DS->Ele[DS->Riv[i].LeftEle-1].BC[j]=-4*(i+1);
				}
			if(DS->Ele[DS->Riv[i].RightEle-1].nabr[j]==DS->Riv[i].LeftEle)
				{
				DS->Ele[DS->Riv[i].RightEle-1].BC[j]=-4*(i+1);
				}
			}
    		DS->Riv[i].x = (DS->Node[DS->Riv[i].FromNode-1].x + DS->Node[DS->Riv[i].ToNode-1].x)/2;
    		DS->Riv[i].y = (DS->Node[DS->Riv[i].FromNode-1].y + DS->Node[DS->Riv[i].ToNode-1].y)/2;
    		DS->Riv[i].zmax = (DS->Node[DS->Riv[i].FromNode-1].zmax + DS->Node[DS->Riv[i].ToNode-1].zmax)/2;
    		DS->Riv[i].depth = CS->Cal.rivDepth*DS->Riv_Shape[DS->Riv[i].shape-1].depth;
		DS->Riv[i].coeff=CS->Cal.rivShapeCoeff*DS->Riv_Shape[DS->Riv[i].shape - 1].coeff;
    		DS->Riv[i].zmin = DS->Riv[i].zmax - DS->Riv[i].depth;  
    		DS->Riv[i].Length = sqrt(pow(DS->Node[DS->Riv[i].FromNode-1].x -DS->Node[DS->Riv[i].ToNode-1].x, 2) + pow(DS->Node[DS->Riv[i].FromNode-1].y - DS->Node[DS->Riv[i].ToNode-1].y, 2));
		DS->Riv[i].KsatH=CS->Cal.rivKsatH*DS->Riv_Mat[DS->Riv[i].material-1].KsatH;
		DS->Riv[i].KsatV=CS->Cal.rivKsatV*DS->Riv_Mat[DS->Riv[i].material-1].KsatV;
		DS->Riv[i].bedThick=CS->Cal.rivbedThick*DS->Riv_Mat[DS->Riv[i].material-1].bedThick;
		DS->Riv[i].Rough=CS->Cal.rivRough*DS->Riv_Mat[DS->Riv[i].material - 1].Rough;
		/* Initialization for rectangular cells beneath river */
		/* Note: Ideally this data should be read from the decomposition itself */
		/* but it is not supported right now in PIHMgis (Bhatt, G and Kumar, M; 2007) */
		DS->Ele[i+DS->NumEle].zmax=DS->Riv[i].zmin;
		DS->Ele[i+DS->NumEle].zmin=DS->Riv[i].zmax-(0.5*(DS->Ele[DS->Riv[i].LeftEle-1].zmax+DS->Ele[DS->Riv[i].RightEle-1].zmax)-0.5*(DS->Ele[DS->Riv[i].LeftEle-1].zmin+DS->Ele[DS->Riv[i].RightEle-1].zmin));
//		DS->Ele[i+DS->NumEle].zmin=DS->Riv[i].zmax-40;
		DS->Ele[i+DS->NumEle].macD=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macD+DS->Ele[DS->Riv[i].RightEle-1].macD)>DS->Riv[i].depth?0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macD+DS->Ele[DS->Riv[i].RightEle-1].macD)-DS->Riv[i].depth:0;	
		DS->Ele[i+DS->NumEle].macKsatH=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macKsatH+DS->Ele[DS->Riv[i].RightEle-1].macKsatH);
		DS->Ele[i+DS->NumEle].vAreaF=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].vAreaF+DS->Ele[DS->Riv[i].RightEle-1].vAreaF);
		DS->Ele[i+DS->NumEle].KsatH=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].KsatH+DS->Ele[DS->Riv[i].RightEle-1].KsatH);
		DS->Ele[i+DS->NumEle].Porosity=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].Porosity+DS->Ele[DS->Riv[i].RightEle-1].Porosity);
 		}
	for(i=0;i<DS->NumPrep;i++)
		{
		for(j=0; j<DS->TSD_Prep[i].length; j++)
			{
			DS->TSD_Prep[i].TS[j][1]=CS->Cal.Prep*DS->TSD_Prep[i].TS[j][1];
			}
		}
        for(i=0;i<DS->NumTemp;i++)
                {
                for(j=0; j<DS->TSD_Temp[i].length; j++)
                        {
                        DS->TSD_Temp[i].TS[j][1]=CS->Cal.Temp*DS->TSD_Temp[i].TS[j][1];
                        }
                }
	/* Memory allocation of print variables */	
        for(i=0;i<22;i++)
                {
                if(i==0)
                        {
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumEle+DS->NumRiv,sizeof(realtype));
                        }
                else if((i>=7)&&(i<19))
                        {
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumRiv,sizeof(realtype));
                        }
                else
                        {
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumEle,sizeof(realtype));
                        }
                }
	/* Debugging artifacts in data created due to coarser resolution of model elements */
	if(CS->Debug==1)
		{
		for(i=0; i<DS->NumEle; i++)
			{
			/* Correction of Surf Elev (artifacts due to coarse scale discretization). Not needed if there is lake feature.*/
			tmpBool=1;
			for(j=0;j<3;j++)
				{
				if(DS->Ele[i].nabr[j]>0)
					{
					tempvalue1=DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax;
					if(DS->Ele[i].zmax-tempvalue1>=0)
						{
						tmpBool=0;
						break;
						}
					}
				}	
			if(tmpBool==1)
				{
				printf("\n Ele %d is sink ",i+1);
				/* Note: Following correction is being applied for debug==1 case only */
				printf("\tBfore: %lf Corrected using:",DS->Ele[i].zmax); 
				tempvalue1=10000000;
				for(j=0;j<3;j++)
					{
					if(DS->Ele[i].nabr[j]>0)
						{
						DS->Ele[i].zmax=(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax);
						tempvalue1=tempvalue1>DS->Ele[i].zmax?DS->Ele[i].zmax:tempvalue1;
						printf("(%d)%lf  ",j+1,(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax));
						}
					}
				DS->Ele[i].zmax=tempvalue1;
				printf("=(New)%lf  ",DS->Ele[i].zmax);
				} 
			}
		/* Correction of BedRck Elev. Is this needed? */ 
		printf("\n Do you want to correct Bed Rock Elev too (1[y]/0[n])");
		scanf("%d",&BoolBR);
		if(BoolBR==1)
		{
		for(i=0; i<DS->NumEle; i++)
			{
			tmpBool=1;
			for(j=0;j<3;j++)
				{
				if(DS->Ele[i].nabr[j]>0)
					{
					tempvalue1=DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmin:DS->Ele[-(DS->Ele[i].BC[j]/4)-1+DS->NumEle].zmin;
					if(DS->Ele[i].zmin-tempvalue1>=0)
						{
						tmpBool=0;
						break;
						}
					}
				}	
			if(tmpBool==1)
				{
				printf("\n Ele %d is sink ",i+1);
				/* Note: Following correction is being applied for debug==1 case only */
				printf("\tBfore: %lf Corrected using:",DS->Ele[i].zmin); 
				tempvalue1=10000000;
				for(j=0;j<3;j++)
					{
					if(DS->Ele[i].nabr[j]>0)
						{
						DS->Ele[i].zmin=(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmin:DS->Ele[-(DS->Ele[i].BC[j]/4)-1+DS->NumEle].zmin);
						tempvalue1=tempvalue1>DS->Ele[i].zmin?DS->Ele[i].zmin:tempvalue1;
						printf("(%d)%lf  ",j+1,(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmin:DS->Ele[-(DS->Ele[i].BC[j]/4)-1+DS->NumEle].zmin));
						}
					}
				DS->Ele[i].zmin=tempvalue1;
				printf("=(New)%lf  ",DS->Ele[i].zmin);
				} 
			}
		}
		getchar();
		printf("\nHit any key to see more details");
		for(i=0;i<DS->NumRiv;i++)
			{
			if(DS->Riv[i].down>0)
				{
				if(DS->Riv[i].zmin<DS->Riv[DS->Riv[i].down-1].zmin)
					{
					BoolR=1;
					printf("\n Riv %d is lower than downstream Riv %d",i+1,DS->Riv[i].down);
					}
				}
			}
		if(BoolR==1)
			{
			printf("\n\tRiver elevation correction needed");
			getchar();
			}	
		}
  	for(i=0; i<DS->NumEle; i++)
  		{
                a_x = DS->Node[DS->Ele[i].node[0]-1].x;
                b_x = DS->Node[DS->Ele[i].node[1]-1].x;
                c_x = DS->Node[DS->Ele[i].node[2]-1].x;
                a_y = DS->Node[DS->Ele[i].node[0]-1].y;
                b_y = DS->Node[DS->Ele[i].node[1]-1].y;
                c_y = DS->Node[DS->Ele[i].node[2]-1].y;
		for(j=0;j<3;j++)
			{
			/* Note: Assumption here is that the forumulation is circumcenter based */
			switch(j)
				{
				case 0:
                			distX=(DS->Ele[i].x-0.5*(b_x+c_x));
                			distY=(DS->Ele[i].y-0.5*(b_y+c_y));
					break;
                                case 1:
                			distX=(DS->Ele[i].x-0.5*(c_x+a_x));
                			distY=(DS->Ele[i].y-0.5*(c_y+a_y));
					break;
                                case 2:
                			distX=(DS->Ele[i].x-0.5*(a_x+b_x));
                			distY=(DS->Ele[i].y-0.5*(a_y+b_y));
					break;
				}
			DS->Ele[i].surfH[j]=(DS->Ele[i].nabr[j]>0)?(DS->Ele[i].BC[j]>-4?(DS->Ele[DS->Ele[i].nabr[j]-1].zmax):DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax):DS->Ele[i].BC[j]<=-4?DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax:(DS->Ele[i].zmax); 
			DS->Ele[i].surfX[j]=(DS->Ele[i].nabr[j]>0)?(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].x:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].x):(DS->Ele[i].x-2*distX);          
			DS->Ele[i].surfY[j]=DS->Ele[i].nabr[j]>0?(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].y:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].y):(DS->Ele[i].y-2*distY);
                	}
		DS->Ele[i].dhBYdx=-(DS->Ele[i].surfY[2]*(DS->Ele[i].surfH[1]-DS->Ele[i].surfH[0])+DS->Ele[i].surfY[1]*(DS->Ele[i].surfH[0]-DS->Ele[i].surfH[2])+DS->Ele[i].surfY[0]*(DS->Ele[i].surfH[2]-DS->Ele[i].surfH[1]))/(DS->Ele[i].surfX[2]*(DS->Ele[i].surfY[1]-DS->Ele[i].surfY[0])+DS->Ele[i].surfX[1]*(DS->Ele[i].surfY[0]-DS->Ele[i].surfY[2])+DS->Ele[i].surfX[0]*(DS->Ele[i].surfY[2]-DS->Ele[i].surfY[1]));  
		DS->Ele[i].dhBYdy=-(DS->Ele[i].surfX[2]*(DS->Ele[i].surfH[1]-DS->Ele[i].surfH[0])+DS->Ele[i].surfX[1]*(DS->Ele[i].surfH[0]-DS->Ele[i].surfH[2])+DS->Ele[i].surfX[0]*(DS->Ele[i].surfH[2]-DS->Ele[i].surfH[1]))/(DS->Ele[i].surfY[2]*(DS->Ele[i].surfX[1]-DS->Ele[i].surfX[0])+DS->Ele[i].surfY[1]*(DS->Ele[i].surfX[0]-DS->Ele[i].surfX[2])+DS->Ele[i].surfY[0]*(DS->Ele[i].surfX[2]-DS->Ele[i].surfX[1]));
  		}
  	/* initialize state variable */
  	/* relax case */
	if (CS->init_type == 0)
  		{
    		for(i=0; i<DS->NumEle; i++)
    			{
      			DS->EleIS[i] = 0;
			DS->EleSnow[i]=0;
			/* Note Two components can be separately read too */
			DS->EleSnowGrnd[i]=(1-DS->Ele[i].VegFrac)* DS->EleSnow[i];
			DS->EleSnowCanopy[i]=DS->Ele[i].VegFrac* DS->EleSnow[i];
      			NV_Ith_S(CV_Y, i) = 0;
      			NV_Ith_S(CV_Y, i + DS->NumEle) = 0;
      			NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele[i].zmax - DS->Ele[i].zmin -0.1;
			}
    		for(i=0; i<DS->NumRiv; i++)
    			{
      			NV_Ith_S(CV_Y, i + 3*DS->NumEle) = 0;
			/* Note once the element beneath river is incorporated in decomposition and .mesh file, initialization should be perfomed based on the location data instead of average of neighbor properties */
			NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv)=(DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin) -0.1;
    			}
  		}
  	/* data initialization mode */  
  	else if (CS->init_type == 1)
  		{
		if(DS->UnsatMode ==1)
			{    
    			}
		if(DS->UnsatMode ==2)
			{    
    			for(i=0; i<DS->NumEle; i++)
    				{
      				DS->EleIS[i] = DS->Ele_IC[i].interception;
				DS->EleSnow[i]=DS->Ele_IC[i].snow;
				/* Note Two components can be separately read too */
				DS->EleSnowGrnd[i]=(1-DS->Ele[i].VegFrac)* DS->EleSnow[i];
				DS->EleSnowCanopy[i]=DS->Ele[i].VegFrac* DS->EleSnow[i];
      				NV_Ith_S(CV_Y, i) = DS->Ele_IC[i].surf;
				/* Note: delete 0.1 here */
      				NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[i].unsat;
				NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele_IC[i].sat; 
				/* Note: delete line below for general*/
//     				NV_Ith_S(CV_Y, i + 2*DS->NumEle) = 0*DS->Ele_IC[i].sat+(DS->Ele[i].zmax - DS->Ele[i].zmin)*0.1;
      				if ((NV_Ith_S(CV_Y, i + DS->NumEle) + NV_Ith_S(CV_Y, i + 2*DS->NumEle)) >= (DS->Ele[i].zmax - DS->Ele[i].zmin))
      					{
        				NV_Ith_S(CV_Y, i + DS->NumEle) = ((DS->Ele[i].zmax - DS->Ele[i].zmin) - NV_Ith_S(CV_Y, i + 2*DS->NumEle))*0.98;
        				if (NV_Ith_S(CV_Y, i + DS->NumEle) < 0) 
						{
						NV_Ith_S(CV_Y, i + DS->NumEle) = 0; 
						}
      					} 
    				}	  
    			for(i=0; i<DS->NumRiv; i++)
    				{
      				NV_Ith_S(CV_Y, i + 3*DS->NumEle) = DS->Riv_IC[DS->Riv[i].IC-1].value;
				/* Note once the element beneath river is incorporated in decomposition and .mesh file, initialization should be perfomed based on the location data instead of average of neighbor properties */
//				NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv) = 0.5*(DS->Ele_IC[DS->Riv[i].LeftEle-1].sat+DS->Ele_IC[DS->Riv[i].RightEle-1].sat);
				NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv) = (DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin)-0.1;
    				}
    			}    	
  		}  
  	/* hot start mode */
  	else
  		{
    		fn = (char *)malloc((strlen(filename)+5)*sizeof(char));
    		strcpy(fn, filename);
    		init_file = fopen(strcat(fn, ".init"), "r");
  
    		if(init_file == NULL)
    			{
      			printf("\n  Fatal Error: %s.init is in use or does not exist!\n", filename);
      			exit(1);
    			}
    		else
    			{
      			for(i=0; i<DS->NumEle; i++)
      				{
        			fscanf(init_file, "%lf %lf %lf %lf %lf", &DS->EleIS[i],&DS->EleSnow[i],&tempvalue1,&tempvalue2,&tempvalue3);
                                DS->EleSnowGrnd[i]=(1-DS->Ele[i].VegFrac)* DS->EleSnow[i];
                                DS->EleSnowCanopy[i]=DS->Ele[i].VegFrac* DS->EleSnow[i];
				NV_Ith_S(CV_Y, i)=tempvalue1;
				NV_Ith_S(CV_Y, i + DS->NumEle)=tempvalue2;
				NV_Ith_S(CV_Y, i + 2*DS->NumEle)=tempvalue3;	
      				}  
      			for(i=0; i<DS->NumRiv; i++)
      				{
        			fscanf(init_file, "%lf %lf",&tempvalue1,&tempvalue2);
        			NV_Ith_S(CV_Y, i + 3*DS->NumEle) = tempvalue1;
				NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv) =tempvalue2;
      				} 
    			}
    		fclose(init_file); 
  		}
  	printf("done.\n");
	}

