/*********************************************************************************
 * File        : read_alloc.c                                                    *
 * Function    : read parameter files for PIHM 2.0                     	         *
 * Version     : Nov, 2007 (2.0)                                                 *
 * Developer of PIHM2.0:	Mukesh Kumar (muk139@psu.edu)		         * 
 * Developer of PIHM1.0:	Yizhong Qu   (quyizhong@gmail.com)	         * 
 *-------------------------------------------------------------------------------*
 *                                                                               *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0..............................*
 * a) Addition of three new input files: file.calib, file.lc and file.geol       *
 * b) Declaration and allocation  of new variables for new process, shape 	 *
 *    representations  and calibration (e.g. ET, Infiltration, Macropore, 	 *
 *    Stormflow, Element beneath river, river shapes, river bed property, 	 *
 *    thresholds for root zone, infiltration and macropore depths, land cover    * 
 *    attributes etc)                                                            *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                         *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sundialstypes.h"
#include "pihm.h"  


void read_alloc(char *filename, Model_Data DS, Control_Data *CS)
	{
  	int i, j;
  	int tempindex;
  
  	int NumTout;
  	char *fn[9];
  	char tempchar[5];
  
  	FILE *mesh_file;	/* Pointer to .mesh file */
  	FILE *att_file;		/* Pointer to .att file */
  	FILE *forc_file;	/* Pointer to .forc file*/
  	FILE *ibc_file;		/* Pointer to .ibc file*/
  	FILE *soil_file;	/* Pointer to .soil file */
  	FILE *geol_file;	/* Pointer to .geol file */
  	FILE *lc_file;		/* Pointer to .lc file */
  	FILE *para_file;	/* Pointer to .para file*/
  	FILE *riv_file;		/* Pointer to .riv file */
	FILE *global_calib;	/* Pointer to .calib file */
 
  
  	printf("\nStart reading in input files ... \n");
  
  	/*========== open *.riv file ==========*/ 
  	printf("\n  1) reading %s.riv  ... ", filename);
  	fn[0] = (char *)malloc((strlen(filename)+4)*sizeof(char));
  	strcpy(fn[0], filename);
  	riv_file =  fopen(strcat(fn[0], ".riv"), "r");
  
  	if(riv_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.riv is in use or does not exist!\n", filename);
    		exit(1);
  		}
  
  	/* start reading riv_file */ 
  	fscanf(riv_file, "%d", &DS->NumRiv);
  
  	DS->Riv = (river_segment *)malloc(DS->NumRiv*sizeof(river_segment));
  	DS->Riv_IC = (river_IC *)malloc(DS->NumRiv*sizeof(river_IC));
  
  	for (i=0; i<DS->NumRiv; i++)
  		{
    		fscanf(riv_file, "%d", &(DS->Riv[i].index));
    		fscanf(riv_file, "%d %d", &(DS->Riv[i].FromNode), &(DS->Riv[i].ToNode));
    		fscanf(riv_file, "%d", &(DS->Riv[i].down));
    		fscanf(riv_file, "%d %d", &(DS->Riv[i].LeftEle), &(DS->Riv[i].RightEle));
    		fscanf(riv_file, "%d %d", &(DS->Riv[i].shape), &(DS->Riv[i].material));
    		fscanf(riv_file, "%d %d", &(DS->Riv[i].IC), &(DS->Riv[i].BC));  
    		fscanf(riv_file, "%d", &(DS->Riv[i].reservoir));                        
  		} 
  
  	fscanf(riv_file, "%s %d", tempchar, &DS->NumRivShape);
  	DS->Riv_Shape = (river_shape *)malloc(DS->NumRivShape*sizeof(river_shape));
  
  	for (i=0; i<DS->NumRivShape; i++)
  		{
    		fscanf(riv_file, "%d", &DS->Riv_Shape[i].index);
    		fscanf(riv_file, "%lf", &DS->Riv_Shape[i].depth);
    		fscanf(riv_file, "%d %lf",&DS->Riv_Shape[i].interpOrd,&DS->Riv_Shape[i].coeff);
  		}
  
  	fscanf(riv_file, "%s %d", tempchar, &DS->NumRivMaterial);
  	DS->Riv_Mat = (river_material *)malloc(DS->NumRivMaterial*sizeof(river_material));
  
  	for (i=0; i<DS->NumRivMaterial; i++)
  		{
    		fscanf(riv_file, "%d %lf %lf %lf %lf %lf", &DS->Riv_Mat[i].index, &DS->Riv_Mat[i].Rough, &DS->Riv_Mat[i].Cwr, &DS->Riv_Mat[i].KsatH,&DS->Riv_Mat[i].KsatV,&DS->Riv_Mat[i].bedThick);
  		}
  
  	fscanf(riv_file, "%s %d", tempchar, &DS->NumRivIC);
  	DS->Riv_IC = (river_IC *)malloc(DS->NumRivIC*sizeof(river_IC));
  
  	for (i=0; i<DS->NumRivIC; i++)
  		{
    		fscanf(riv_file, "%d %lf", &DS->Riv_IC[i].index, &DS->Riv_IC[i].value);
  		}

  	fscanf(riv_file, "%s %d", tempchar, &DS->NumRivBC);
  	DS->TSD_Riv = (TSD *)malloc(DS->NumRivBC*sizeof(TSD));
  
  	for(i=0; i<DS->NumRivBC; i++)
  		{
    		fscanf(riv_file, "%s %d %d", DS->TSD_Riv[i].name, &DS->TSD_Riv[i].index, &DS->TSD_Riv[i].length);
    
    		DS->TSD_Riv[i].TS = (realtype **)malloc((DS->TSD_Riv[i].length)*sizeof(realtype));
    		for(j=0; j<DS->TSD_Riv[i].length; j++)
    			{
      			DS->TSD_Riv[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Riv[i].length; j++)
    			{
      			fscanf(riv_file, "%lf %lf", &DS->TSD_Riv[i].TS[j][0], &DS->TSD_Riv[i].TS[j][1]);
    			}
  		}  
  
  	// read in reservoir information
  	fscanf(riv_file, "%s %d", tempchar, &DS->NumRes);
  	if(DS->NumRes > 0)
  		{
    		/* read in reservoir information */
    
  		}
  
  	fclose(riv_file);
  	printf("done.\n");
	
  	/*========== open *.mesh file ==========*/
  	printf("\n  2) reading %s.mesh ... ", filename);
  	fn[1] = (char *)malloc((strlen(filename)+5)*sizeof(char));
  	strcpy(fn[1], filename);
  	mesh_file = fopen(strcat(fn[1], ".mesh"), "r");

  	if(mesh_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
    		exit(1);
  		}
    
  	/* start reading mesh_file */ 
  	fscanf(mesh_file,"%d %d", &DS->NumEle, &DS->NumNode);
  
  	DS->Ele = (element *)malloc((DS->NumEle+DS->NumRiv)*sizeof(element));
  	DS->Node = (nodes *)malloc(DS->NumNode*sizeof(nodes));
  
  	/* read in elements information */ 
  	for (i=0; i<DS->NumEle; i++)
  		{
    		fscanf(mesh_file, "%d", &(DS->Ele[i].index));
    		fscanf(mesh_file, "%d %d %d", &(DS->Ele[i].node[0]), &(DS->Ele[i].node[1]), &(DS->Ele[i].node[2]));
    		fscanf(mesh_file, "%d %d %d", &(DS->Ele[i].nabr[0]), &(DS->Ele[i].nabr[1]), &(DS->Ele[i].nabr[2]));
  		}
  
  	/* read in nodes information */   
  	for (i=0; i<DS->NumNode; i++)
  		{
    		fscanf(mesh_file, "%d", &(DS->Node[i].index));
    		fscanf(mesh_file, "%lf %lf", &(DS->Node[i].x), &(DS->Node[i].y));
    		fscanf(mesh_file, "%lf %lf", &(DS->Node[i].zmin),&(DS->Node[i].zmax));
  		}  
  
  	printf("done.\n");
  
  	/* finish reading mesh_files */  
  	fclose(mesh_file);
  
  	/*========== open *.att file ==========*/
  	printf("\n  3) reading %s.att  ... ", filename);
  	fn[2] = (char *)malloc((strlen(filename)+4)*sizeof(char));
  	strcpy(fn[2], filename);
  	att_file = fopen(strcat(fn[2], ".att"), "r");

  	if(att_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.att is in use or does not exist!\n", filename);
    		exit(1);
  		}
    
  	/* start reading att_file */ 
  	DS->Ele_IC = (element_IC *)malloc(DS->NumEle*sizeof(element_IC));
  	for (i=0; i<DS->NumEle; i++)
  		{
    		fscanf(att_file, "%d", &(tempindex));
    		fscanf(att_file, "%d %d %d", &(DS->Ele[i].soil), &(DS->Ele[i].geol), &(DS->Ele[i].LC));
    		fscanf(att_file, "%lf %lf %lf %lf %lf",&(DS->Ele_IC[i].interception),&(DS->Ele_IC[i].snow),&(DS->Ele_IC[i].surf),&(DS->Ele_IC[i].unsat),&(DS->Ele_IC[i].sat));
    		fscanf(att_file, "%d %d", &(DS->Ele[i].prep), &(DS->Ele[i].temp));
    		fscanf(att_file, "%d %d", &(DS->Ele[i].humidity), &(DS->Ele[i].WindVel));
    		fscanf(att_file, "%d %d", &(DS->Ele[i].Rn), &(DS->Ele[i].G));
    		fscanf(att_file, "%d %d %d", &(DS->Ele[i].pressure), &(DS->Ele[i].source), &(DS->Ele[i].meltF));
		for(j=0;j<3;j++)
			{
    			fscanf(att_file, "%d", &(DS->Ele[i].BC[j]));
			}
		fscanf(att_file, "%d", &(DS->Ele[i].Macropore));
 		}
  
  	printf("done.\n");
  
  	/* finish reading att_files */  
  	fclose(att_file);
    
  	/*========== open *.soil file ==========*/  
  	printf("\n  4) reading %s.soil ... ", filename);
  	fn[3] = (char *)malloc((strlen(filename)+5)*sizeof(char));
  	strcpy(fn[3], filename);
  	soil_file = fopen(strcat(fn[3], ".soil"), "r");
  
  	if(soil_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.soil is in use or does not exist!\n", filename);
    		exit(1);
  		}
  
  	/* start reading soil_file */  
  	fscanf(soil_file, "%d", &DS->NumSoil);
  	DS->Soil = (soils *)malloc(DS->NumSoil*sizeof(soils));
  
  	for (i=0; i<DS->NumSoil; i++)
  		{
    		fscanf(soil_file, "%d", &(DS->Soil[i].index));
		/* Note: Soil KsatH and macKsatH is not used in model calculation anywhere */
    		fscanf(soil_file, "%lf",&(DS->Soil[i].KsatV));
    		fscanf(soil_file, "%lf %lf %lf", &(DS->Soil[i].ThetaS), &(DS->Soil[i].ThetaR), &(DS->Soil[i].infD));
    		fscanf(soil_file, "%lf %lf", &(DS->Soil[i].Alpha), &(DS->Soil[i].Beta));
    		fscanf(soil_file, "%lf %lf", &(DS->Soil[i].hAreaF),&(DS->Soil[i].macKsatV));
  		} 
 
  	fclose(soil_file);
  	printf("done.\n");
      
        /*========== open *.geol file ==========*/
        printf("\n  5) reading %s.geol ... ", filename);
        fn[4] = (char *)malloc((strlen(filename)+5)*sizeof(char));
        strcpy(fn[4], filename);
        geol_file = fopen(strcat(fn[4], ".geol"), "r");

        if(geol_file == NULL)
                {
                printf("\n  Fatal Error: %s.geol is in use or does not exist!\n", filename);
                exit(1);
                }

        /* start reading soil_file */
        fscanf(geol_file, "%d", &DS->NumGeol);
        DS->Geol = (geol *)malloc(DS->NumGeol*sizeof(geol));

        for (i=0; i<DS->NumGeol; i++)
                {
                fscanf(geol_file, "%d", &(DS->Geol[i].index));
                /* Geol macKsatV is not used in model calculation anywhere */
                fscanf(geol_file, "%lf %lf", &(DS->Geol[i].KsatH),&(DS->Geol[i].KsatV));
                fscanf(geol_file, "%lf %lf", &(DS->Geol[i].ThetaS), &(DS->Geol[i].ThetaR));
                fscanf(soil_file, "%lf %lf", &(DS->Geol[i].Alpha), &(DS->Geol[i].Beta));
                fscanf(soil_file, "%lf %lf %lf", &(DS->Geol[i].vAreaF),&(DS->Geol[i].macKsatH),&(DS->Geol[i].macD));
                }

        fclose(geol_file);
        printf("done.\n");


  	/*========== open *.lc file ==========*/  
  	printf("\n  6) reading %s.lc ... ", filename);
  	fn[5] = (char *)malloc((strlen(filename)+3)*sizeof(char));
  	strcpy(fn[5], filename);
  	lc_file = fopen(strcat(fn[5], ".lc"), "r");
  
  	if(lc_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.land cover is in use or does not exist!\n", filename);
    		exit(1);
  		}
  
  	/* start reading land cover file */  
  	fscanf(lc_file, "%d", &DS->NumLC);
  
  	DS->LandC = (LC *)malloc(DS->NumLC*sizeof(LC));
  
  	for (i=0; i<DS->NumLC; i++)
  		{
    		fscanf(lc_file, "%d", &(DS->LandC[i].index));
    		fscanf(lc_file, "%lf", &(DS->LandC[i].LAImax));
    		fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Rmin), &(DS->LandC[i].Rs_ref));
    		fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Albedo), &(DS->LandC[i].VegFrac));
    		fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Rough),&(DS->LandC[i].RzD));
  		} 
 
  	fclose(lc_file);
  	printf("done.\n");

  
  	/*========== open *.forc file ==========*/  
  	printf("\n  7) reading %s.forc ... ", filename);
  	fn[6] = (char *)malloc((strlen(filename)+5)*sizeof(char));
  	strcpy(fn[6], filename);
  	forc_file = fopen(strcat(fn[6], ".forc"), "r");
  
  	if(forc_file == NULL)
  	{
    	printf("\n  Fatal Error: %s.forc is in use or does not exist!\n", filename);
    	exit(1);
  	}
  
  	/* start reading forc_file */
  	fscanf(forc_file, "%d %d", &DS->NumPrep, &DS->NumTemp);
  	fscanf(forc_file, "%d %d", &DS->NumHumidity, &DS->NumWindVel);
  	fscanf(forc_file, "%d %d", &DS->NumRn, &DS->NumG);
  	fscanf(forc_file, "%d %d", &DS->NumP, &DS->NumLC);
  	fscanf(forc_file, "%d", &DS->NumMeltF);
  	fscanf(forc_file, "%d", &DS->NumSource);
  
  	DS->TSD_Prep = (TSD *)malloc(DS->NumPrep*sizeof(TSD));
  	DS->TSD_Temp = (TSD *)malloc(DS->NumTemp*sizeof(TSD));
  	DS->TSD_Humidity = (TSD *)malloc(DS->NumHumidity*sizeof(TSD));
  	DS->TSD_WindVel = (TSD *)malloc(DS->NumWindVel*sizeof(TSD));
  	DS->TSD_Rn = (TSD *)malloc(DS->NumRn*sizeof(TSD));
  	DS->TSD_G = (TSD *)malloc(DS->NumG*sizeof(TSD));
  	DS->TSD_Pressure = (TSD *)malloc(DS->NumP*sizeof(TSD));
  	DS->TSD_LAI = (TSD *)malloc(DS->NumLC*sizeof(TSD));
  	DS->TSD_RL = (TSD *)malloc(DS->NumLC*sizeof(TSD));  
  	DS->TSD_MeltF = (TSD *)malloc(DS->NumMeltF*sizeof(TSD));
  	DS->TSD_Source = (TSD *)malloc(DS->NumSource*sizeof(TSD));
  
  	DS->ISFactor = (realtype *)malloc(DS->NumLC*sizeof(realtype));
  	DS->windH = (realtype *)malloc(DS->NumWindVel*sizeof(realtype));
     
  	for(i=0; i<DS->NumPrep; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_Prep[i].name, &DS->TSD_Prep[i].index, &DS->TSD_Prep[i].length);
    
    		DS->TSD_Prep[i].TS = (realtype **)malloc((DS->TSD_Prep[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_Prep[i].length; j++)
    			{
      			DS->TSD_Prep[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Prep[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_Prep[i].TS[j][0], &DS->TSD_Prep[i].TS[j][1]);
    			}
  		}  
  
  	for(i=0; i<DS->NumTemp; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_Temp[i].name, &DS->TSD_Temp[i].index, &DS->TSD_Temp[i].length);
    
    		DS->TSD_Temp[i].TS = (realtype **)malloc((DS->TSD_Temp[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_Temp[i].length; j++)
    			{
      			DS->TSD_Temp[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Temp[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_Temp[i].TS[j][0], &DS->TSD_Temp[i].TS[j][1]);
    			}
  		} 
  
  	for(i=0; i<DS->NumHumidity; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_Humidity[i].name, &DS->TSD_Humidity[i].index, &DS->TSD_Humidity[i].length);
    
    		DS->TSD_Humidity[i].TS = (realtype **)malloc((DS->TSD_Humidity[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_Humidity[i].length; j++)
    			{
      			DS->TSD_Humidity[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Humidity[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_Humidity[i].TS[j][0], &DS->TSD_Humidity[i].TS[j][1]);
    			}
  		} 
  
  	for(i=0; i<DS->NumWindVel; i++)
  		{
    		fscanf(forc_file, "%s %d %d %lf", DS->TSD_WindVel[i].name, &DS->TSD_WindVel[i].index, &DS->TSD_WindVel[i].length, &DS->windH[i]);
    		DS->TSD_WindVel[i].TS = (realtype **)malloc((DS->TSD_WindVel[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_WindVel[i].length; j++)
    			{
      			DS->TSD_WindVel[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_WindVel[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_WindVel[i].TS[j][0], &DS->TSD_WindVel[i].TS[j][1]);
    			}
  		} 

  	for(i=0; i<DS->NumRn; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_Rn[i].name, &DS->TSD_Rn[i].index, &DS->TSD_Rn[i].length);
    
    		DS->TSD_Rn[i].TS = (realtype **)malloc((DS->TSD_Rn[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_Rn[i].length; j++)
    			{
      			DS->TSD_Rn[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Rn[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_Rn[i].TS[j][0], &DS->TSD_Rn[i].TS[j][1]);
    			}
  		} 

  	for(i=0; i<DS->NumG; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_G[i].name, &DS->TSD_G[i].index, &DS->TSD_G[i].length);
    
    		DS->TSD_G[i].TS = (realtype **)malloc((DS->TSD_G[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_G[i].length; j++)
    			{
      			DS->TSD_G[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_G[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_G[i].TS[j][0], &DS->TSD_G[i].TS[j][1]);
    			}
  		} 

  	for(i=0; i<DS->NumP; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_Pressure[i].name, &DS->TSD_Pressure[i].index, &DS->TSD_Pressure[i].length);
    
    		DS->TSD_Pressure[i].TS = (realtype **)malloc((DS->TSD_Pressure[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_Pressure[i].length; j++)
    			{
      			DS->TSD_Pressure[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Pressure[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_Pressure[i].TS[j][0], &DS->TSD_Pressure[i].TS[j][1]);
    			}
  		} 

  	for(i=0; i<DS->NumLC; i++)
  		{
    		fscanf(forc_file, "%s %d %d %lf", DS->TSD_LAI[i].name, &DS->TSD_LAI[i].index, &DS->TSD_LAI[i].length, &DS->ISFactor[i]);
    
    		DS->TSD_LAI[i].TS = (realtype **)malloc((DS->TSD_LAI[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_LAI[i].length; j++)
    			{
      			DS->TSD_LAI[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_LAI[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_LAI[i].TS[j][0], &DS->TSD_LAI[i].TS[j][1]);
    			}
  		} 

  	for(i=0; i<DS->NumLC; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_RL[i].name, &DS->TSD_RL[i].index, &DS->TSD_RL[i].length);
    
    		DS->TSD_RL[i].TS = (realtype **)malloc((DS->TSD_RL[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_RL[i].length; j++)
    			{
      			DS->TSD_RL[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_RL[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_RL[i].TS[j][0], &DS->TSD_RL[i].TS[j][1]);
    			}
  		} 
  
  	for(i=0; i<DS->NumMeltF; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_MeltF[i].name, &DS->TSD_MeltF[i].index, &DS->TSD_MeltF[i].length);
    
    		DS->TSD_MeltF[i].TS = (realtype **)malloc((DS->TSD_MeltF[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_MeltF[i].length; j++)
    			{
      			DS->TSD_MeltF[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_MeltF[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_MeltF[i].TS[j][0], &DS->TSD_MeltF[i].TS[j][1]);
    			}
  		} 
  
  	for(i=0; i<DS->NumSource; i++)
  		{
    		fscanf(forc_file, "%s %d %d", DS->TSD_Source[i].name, &DS->TSD_Source[i].index, &DS->TSD_Source[i].length);
    
    		DS->TSD_Source[i].TS = (realtype **)malloc((DS->TSD_Source[i].length)*sizeof(realtype));
    
    		for(j=0; j<DS->TSD_Source[i].length; j++)
    			{
      			DS->TSD_Source[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
    			}
    
    		for(j=0; j<DS->TSD_Source[i].length; j++)
    			{
      			fscanf(forc_file, "%lf %lf", &DS->TSD_Source[i].TS[j][0], &DS->TSD_Source[i].TS[j][1]);
    			}
  		} 

  	fclose(forc_file);
  	printf("done.\n");
      
  	/*========== open *.ibc file ==========*/     
  	printf("\n  8) reading %s.ibc  ... ", filename);  
  	fn[7] = (char *)malloc((strlen(filename)+4)*sizeof(char));
  	strcpy(fn[7], filename);
  	ibc_file =  fopen(strcat(fn[7], ".ibc"), "r");
  
  	if(ibc_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.ibc is in use or does not exist!\n", filename);
    		exit(1);
  		}
  
  	/* start reading ibc_file */
  	fscanf(ibc_file, "%d %d", &DS->Num1BC, &DS->Num2BC);
  
  	if(DS->Num1BC+DS->Num2BC > 0)
  		{
    		DS->TSD_EleBC = (TSD *)malloc((DS->Num1BC+DS->Num2BC)*sizeof(TSD));
  		}
  
  	if(DS->Num1BC>0)
  		{
    		/* For elements with Dirichilet Boundary Conditions */
    		for(i=0; i<DS->Num1BC; i++)
    			{
      			fscanf(ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index,&DS->TSD_EleBC[i].length);
      
      			DS->TSD_EleBC[i].TS = (realtype **)malloc((DS->TSD_EleBC[i].length)*sizeof(realtype));
      
      			for(j=0; j<DS->TSD_EleBC[i].length; j++)
      				{
        			DS->TSD_EleBC[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
      				}
      
      			for(j=0; j<DS->TSD_EleBC[i].length; j++)
      				{
        			fscanf(forc_file, "%lf %lf", &DS->TSD_EleBC[i].TS[j][0],&DS->TSD_EleBC[i].TS[j][1]);
      				}
    			}    
  		}
  
  	if(DS->Num2BC>0)
  		{
    		/* For elements with Neumann (non-natural) Boundary Conditions */
    		for(i=DS->Num1BC; i<DS->Num1BC+DS->Num2BC; i++)
    			{
      			fscanf(ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index,&DS->TSD_EleBC[i].length);
      
      			DS->TSD_EleBC[i].TS = (realtype **)malloc((DS->TSD_EleBC[i].length)*sizeof(realtype));
      
      			for(j=0; j<DS->TSD_EleBC[i].length; j++)
      				{
        			DS->TSD_EleBC[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
      				}
      			for(j=0; j<DS->TSD_EleBC[i].length; j++)
      				{
        			fscanf(forc_file, "%lf %lf", &DS->TSD_EleBC[i].TS[j][0],&DS->TSD_EleBC[i].TS[j][1]);
      				}
    			}     
  		}
  	fclose(ibc_file);
  	printf("done.\n");
  
  	/*========== open *.para file ==========*/ 
  	printf("\n  9) reading %s.para ... ", filename); 
  	fn[8] = (char *)malloc((strlen(filename)+5)*sizeof(char));
  	strcpy(fn[8], filename);
  	para_file = fopen(strcat(fn[8], ".para"), "r");  
  
  	if(para_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.para is in use or does not exist!\n", filename);
    		exit(1);
  		}
  
  	/* start reading para_file */
  	fscanf(para_file, "%d %d", &(CS->Verbose), &(CS->Debug));
  	fscanf(para_file, "%d", &(CS->init_type));
  	fscanf(para_file, "%d %d %d %d", &(CS->gwD), &(CS->surfD), &(CS->snowD), &(CS->rivStg));
  	fscanf(para_file, "%d %d %d", &(CS->Rech), &(CS->IsD), &(CS->usD));
	fscanf(para_file, "%d %d %d", &(CS->et[0]), &(CS->et[1]), &(CS->et[2]));
	for(i=0;i<10;i++)
		{
		fscanf(para_file, "%d", &(CS->rivFlx[i]));
		}
        fscanf(para_file, "%d %d %d %d", &(CS->gwDInt), &(CS->surfDInt), &(CS->snowDInt), &(CS->rivStgInt));
        fscanf(para_file, "%d %d %d %d", &(CS->RechInt), &(CS->IsDInt), &(CS->usDInt), &(CS->etInt));
	fscanf(para_file, "%d",&(CS->rivFlxInt));
	
  	fscanf(para_file, "%d %d %d", &DS->UnsatMode, &DS->SurfMode, &DS->RivMode);
  	fscanf(para_file, "%d", &(CS->Solver));
  	if(CS->Solver == 2)
  		{
    		fscanf(para_file, "%d %d %lf", &CS->GSType, &CS->MaxK, &CS->delt);
  		}
  	fscanf(para_file, "%lf %lf", &(CS->abstol), &(CS->reltol));
  	fscanf(para_file, "%lf %lf %lf", &(CS->InitStep), &(CS->MaxStep), &(CS->ETStep));
  	fscanf(para_file, "%lf %lf %d", &(CS->StartTime), &(CS->EndTime), &(CS->outtype));
  	if(CS->outtype == 0)
  		{
    		fscanf(para_file, "%lf %lf", &CS->a, &CS->b);
  		}
  
  	if(CS->a != 1.0)
  		{
    		NumTout = (int)(log(1 - (CS->EndTime - CS->StartTime)*(1 -  CS->a)/CS->b)/log(CS->a));
  		}
  	else
  		{
    		if((CS->EndTime - CS->StartTime)/CS->b - ((int) (CS->EndTime - CS->StartTime)/CS->b) > 0)
    			{
      			NumTout = (int) ((CS->EndTime - CS->StartTime)/CS->b);
    			}
    		else
    			{
      			NumTout = (int) ((CS->EndTime - CS->StartTime)/CS->b - 1);
    			}  
  		}
  
  	CS->NumSteps = NumTout + 1;
  
  	CS->Tout = (realtype *)malloc((CS->NumSteps + 1)*sizeof(realtype));
    
  	for(i=0; i<CS->NumSteps+1; i++)
  		{
    		if(i == 0)
    			{
      			CS->Tout[i] = CS->StartTime;
    			}
    		else
    			{
      			CS->Tout[i] = CS->Tout[i-1] + pow(CS->a, i)*CS->b;
    			}  
  		}
  
  	if(CS->Tout[CS->NumSteps] < CS->EndTime)
  		{
    		CS->Tout[CS->NumSteps] = CS->EndTime;
  		}
  
  	fclose(para_file); 
  	printf("done.\n"); 

	printf("\nStart reading in calibration file...\n");

	/*========= open *.calib file ==========*/
	printf("\n  10) reading %s.calib ... ", filename);
        fn[9] = (char *)malloc((strlen(filename)+6)*sizeof(char));
        strcpy(fn[9], filename);
        global_calib = fopen(strcat(fn[9], ".calib"), "r");

        if(global_calib == NULL)
                {
                printf("\n  Fatal Error: %s.calib is in use or does not exist!\n", filename);
                exit(1);
                }

	/* start reading calib_file */
	fscanf(global_calib,"%lf %lf %lf %lf %lf",&CS->Cal.KsatH,&CS->Cal.KsatV,&CS->Cal.infKsatV,&CS->Cal.macKsatH,&CS->Cal.macKsatV);
	fscanf(global_calib,"%lf %lf %lf",&CS->Cal.infD,&CS->Cal.RzD,&CS->Cal.macD);
	fscanf(global_calib,"%lf %lf %lf",&CS->Cal.Porosity,&CS->Cal.Alpha,&CS->Cal.Beta);
	fscanf(global_calib,"%lf %lf",&CS->Cal.vAreaF,&CS->Cal.hAreaF);
	fscanf(global_calib,"%lf %lf %lf",&CS->Cal.VegFrac,&CS->Cal.Albedo,&CS->Cal.Rough);
	fscanf(global_calib,"%lf %lf",&CS->Cal.Prep,&CS->Cal.Temp);
	fscanf(global_calib,"%lf %lf %lf",&DS->pcCal.Et0,&DS->pcCal.Et1,&DS->pcCal.Et2);
	fscanf(global_calib,"%lf %lf %lf %lf",&CS->Cal.rivRough,&CS->Cal.rivKsatH,&CS->Cal.rivKsatV,&CS->Cal.rivbedThick);
	fscanf(global_calib,"%lf %lf",&CS->Cal.rivDepth,&CS->Cal.rivShapeCoeff);
  	printf("done.\n");
  
  	/* finish reading calib file */  
  	fclose(global_calib);
	}

