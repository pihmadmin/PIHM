/*******************************************************************************
 * File        : read_alloc.c                                                  *
 * Function    : read Input Files and allocate memorey for PIHM                *
 * Programmers : Yizhong Qu   @ Pennsylvania State Univeristy                  *
 *               Mukesh Kumar @ Pennsylvania State Univeristy                  *
 *               Gopal Bhatt  @ Pennsylvania State Univeristy                  *
 * Version     : 2.0 (July 10, 2007)                                           *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * In this file, 7 input files are read in and dynamaticlly allocate memory    *
 * based on control parameters in the input files. Please refer to format.pdf  *
 * for details about the format of all 7 input files. All files has the same   *
 * file prefix provided in the command line when starting PIHM. For example,   *
 * "example" will incur reading in example.mesh, example.soil, etc.            *
 *                                                                             *
 * This code is free for users with research purpose only, if appropriate      *
 * citation is refered. However, there is no warranty in any format for this   *
 * product.                                                                    *
 *                                                                             *
 * For questions or comments, please contact the authors of the reference.     *
 * One who want to use it for other consideration may also contact Dr.Duffy    *
 * at cxd11@psu.edu.                                                           *
 *******************************************************************************/

//! @file read_alloc.c input files are read in and dynamatic allocate memory

/* C Header Files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* SUNDIALS Header Files */
#include "sundials_types.h"
#include "pihm.h"
#include "calib.h"




/***************************************************************
    Function reads all the input files
****************************************************************/
void read_alloc(char *filename, Model_Data DS, Control_Data *CS)
//! Function reads all the input files and initializes members of Model and Control Data Structure
/*! \param filename is Identifier of input files
    \param DS is pointer to model data structure
    \param CS is pointer to control data structure
*/
{
    int i, j;
    int tempindex;

    int NumTout; /*!< Enum value TVal1. */
    char *fn[8];
    char tempchar[5];

    FILE *mesh_file;                        /*    Pointer to .mesh file    */
    FILE *att_file;                         /*    Pointer to .att  file    */
    FILE *forc_file;                        /*    Pointer to .forc file    */
    FILE *ibc_file;                         /*    Pointer to .ibc  file    */
    FILE *soil_file;                        /*    Pointer to .soil file    */
    FILE *lc_file;                          /*    Pointer to .lc     file  */
    FILE *para_file;                        /*    Pointer to .para file    */
    FILE *riv_file;                         /*    Pointer to .riv  file    */

    /* CALIBRATION Parameters    */
	realtype roughEle_CALIB;
	realtype roughRiv_CALIB;
	realtype rivCoeff_CALIB;
	realtype rivDepth_CALIB;
	realtype alpha_CALIB;
	realtype set_MP;
	realtype lai_CALIB;
	realtype vegfrac_CALIB;
    realtype albedo_CALIB;

    /*    Initialize CALIBRATION parameters    */
    roughEle_CALIB=setroughEle_CALIB();
    roughRiv_CALIB=setroughRiv_CALIB();
    rivCoeff_CALIB=setrivCoeff_CALIB();
    rivDepth_CALIB=setrivDepth_CALIB();
    alpha_CALIB=setalpha_CALIB();
    set_MP=setset_MP();
    lai_CALIB=setlai_CALIB();
    vegfrac_CALIB=setvegfrac_CALIB();
    albedo_CALIB=setalbedo_CALIB();


    printf("\nStart reading in input files ... \n");

    /****************************************/
    /*========== open *.mesh file ==========*/
    /****************************************/
    printf("\n  1) reading %s.mesh ... ", filename);
    fn[0] = (char *)malloc((strlen(filename)+5)*sizeof(char));
    strcpy(fn[0], filename);
    mesh_file = fopen(strcat(fn[0], ".mesh"), "r");

    if(mesh_file == NULL)
    {
        printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading mesh_file */
    fscanf(mesh_file,"%d %d", &DS->NumEle, &DS->NumNode);

    DS->Ele = (element *)malloc(DS->NumEle*sizeof(element));          /*    Allocate memory for Number of Elements     */
    DS->Node = (nodes *)malloc(DS->NumNode*sizeof(nodes));            /*    Allocate memory for Number of Nodes        */

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

    fclose(mesh_file);
    /* finish reading mesh_files */

    /***************************************/
    /*========== open *.att file ==========*/
    /***************************************/
    printf("\n  2) reading %s.att  ... ", filename);
    fn[1] = (char *)malloc((strlen(filename)+4)*sizeof(char));
    strcpy(fn[1], filename);
    att_file = fopen(strcat(fn[1], ".att"), "r");

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
        fscanf(att_file, "%d %d", &(DS->Ele[i].soil), &(DS->Ele[i].LC));            /*    Read Soil and Land Cover Type               */
        fscanf(att_file, "%lf %lf %lf %lf %lf",&(DS->Ele_IC[i].interception),&(DS->Ele_IC[i].snow),&(DS->Ele_IC[i].surf),&(DS->Ele_IC[i].unsat),&(DS->Ele_IC[i].sat));
                                                                                    /*    Read Initial States                         */
        fscanf(att_file, "%d", &(DS->Ele[i].BC));                                   /*    Read Boundary Condition                     */
        //DS->Ele[i].BC=0;
        fscanf(att_file, "%d %d", &(DS->Ele[i].prep), &(DS->Ele[i].temp));          /*    Read Precipitation & Temperature Class      */
        fscanf(att_file, "%d %d", &(DS->Ele[i].humidity), &(DS->Ele[i].WindVel));   /*    Read Rel. Humidity & Wind Velocity Class    */
        fscanf(att_file, "%d %d", &(DS->Ele[i].Rn), &(DS->Ele[i].G));               /*    Read Solar Radiation Class    + Dummy       */
        fscanf(att_file, "%d %d", &(DS->Ele[i].pressure), &(DS->Ele[i].source));    /*    Read Pressure Class and Source/Sinks        */
    }

    printf("done.\n");

    fclose(att_file);
    /* finish reading mesh_files */

    /****************************************/
    /*========== open *.soil file ==========*/
    /****************************************/
    printf("\n  3) reading %s.soil ... ", filename);
    fn[2] = (char *)malloc((strlen(filename)+5)*sizeof(char));
    strcpy(fn[2], filename);
    soil_file = fopen(strcat(fn[2], ".soil"), "r");

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
        fscanf(soil_file, "%lf", &(DS->Soil[i].Ksat));                                /*    Read Saturated Hydraulic Conductivity    */
        fscanf(soil_file, "%lf %lf", &(DS->Soil[i].SitaS), &(DS->Soil[i].SitaR));     /*    Read Porosity and Residual Porosity      */
        fscanf(soil_file, "%lf %lf", &(DS->Soil[i].Alpha), &(DS->Soil[i].Beta));      /*    Read Soil Parameters Alpha & Beta        */
        fscanf(soil_file, "%d %lf %lf", &(DS->Soil[i].Macropore), &(DS->Soil[i].base), &(DS->Soil[i].gama));
                                                                                      /*    Read Macropore (0/1) Base & Gamma        */
        fscanf(soil_file, "%lf", &(DS->Soil[i].Sf));                                  /*    Read Soil Friction Slope                 */
        fscanf(soil_file, "%lf", &(DS->Soil[i].RzD));                                 /*    Read Root Zone Depth                     */
        fscanf(soil_file, "%d", &(DS->Soil[i].Inf));

        /*     CALIBRATION    */
        DS->Soil[i].Alpha=alpha_CALIB*DS->Soil[i].Alpha;
        DS->Soil[i].Macropore=set_MP;
    }

    fscanf(soil_file, "%d", &DS->NumInc);

    DS->TSD_Inc = (TSD *)malloc(DS->NumInc*sizeof(TSD));

    for(i=0; i<DS->NumInc; i++)
    {
        fscanf(soil_file, "%s %d %d", DS->TSD_Inc[i].name, &DS->TSD_Inc[i].index,
                                  &DS->TSD_Inc[i].length);

        DS->TSD_Inc[i].TS = (realtype **)malloc(DS->TSD_Inc[i].length*sizeof(realtype));
        for(j=0; j<DS->TSD_Inc[i].length; j++)
        {
            DS->TSD_Inc[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
        }

        for(j=0; j<DS->TSD_Inc[i].length; j++)
        {
            fscanf(soil_file, "%lf %lf", &DS->TSD_Inc[i].TS[j][0],
                                   &DS->TSD_Inc[i].TS[j][1]);
        }
    }

    fclose(soil_file);
    printf("done.\n");
    /* Finish reading soil_file */


    /**************************************/
    /*========== open *.lc file ==========*/
    /**************************************/
    printf("\n  3) reading %s.lc ... ", filename);
    fn[3] = (char *)malloc((strlen(filename)+5)*sizeof(char));
    strcpy(fn[3], filename);
    lc_file = fopen(strcat(fn[3], ".lc"), "r");

    if(lc_file == NULL)
    {
        printf("\n  Fatal Error: %s.land cover is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading land cover file */
    fscanf(lc_file, "%d", &DS->NumLC);                                                    /*    Number of Land Cover Classes             */

    DS->LandC = (LC *)malloc(DS->NumLC*sizeof(LC));

    for (i=0; i<DS->NumLC; i++)
    {
        fscanf(lc_file, "%d", &(DS->LandC[i].index));
        fscanf(lc_file, "%lf", &(DS->LandC[i].LAImax));                                   /*    Read Max LAI                             */
        fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Rmin), &(DS->LandC[i].Rs_ref));         /*    Read Min Stomal Resistance & Reference   */
        fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Albedo), &(DS->LandC[i].VegFrac));      /*    Read Albedo and Vegitation Fraction      */
        fscanf(lc_file, "%lf", &(DS->LandC[i].Rough));                                    /*    Read Manning's Roughnes Coefficient      */

        /*    CALIBRATION    */
        DS->LandC[i].VegFrac=vegfrac_CALIB*DS->LandC[i].VegFrac;
        DS->LandC[i].Albedo=albedo_CALIB*DS->LandC[i].Albedo>1.0?1.0:albedo_CALIB*DS->LandC[i].Albedo;
        DS->LandC[i].Rough=roughEle_CALIB*DS->LandC[i].Rough;
    }

    fclose(lc_file);
    printf("done.\n");
    /* Finish reading land cover file */

    /***************************************/
    /*========== open *.riv file ==========*/
    /***************************************/
    printf("\n  4) reading %s.riv  ... ", filename);
    fn[4] = (char *)malloc((strlen(filename)+4)*sizeof(char));
    strcpy(fn[4], filename);
    riv_file =  fopen(strcat(fn[4], ".riv"), "r");

    if(riv_file == NULL)
    {
        printf("\n  Fatal Error: %s.riv is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading .riv File */
    fscanf(riv_file, "%d", &DS->NumRiv);

    DS->Riv = (river_segment *)malloc(DS->NumRiv*sizeof(river_segment));
    DS->Riv_IC = (river_IC *)malloc(DS->NumRiv*sizeof(river_IC));

    for (i=0; i<DS->NumRiv; i++)
    {
        fscanf(riv_file, "%d", &(DS->Riv[i].index));
        fscanf(riv_file, "%d %d", &(DS->Riv[i].FromNode), &(DS->Riv[i].ToNode));         /*    Read From and To Nodes #              */
        fscanf(riv_file, "%d", &(DS->Riv[i].down));                                      /*    Read Down Segment                     */
        fscanf(riv_file, "%d %d", &(DS->Riv[i].LeftEle), &(DS->Riv[i].RightEle));        /*    Read Left and Right Elements #        */
        fscanf(riv_file, "%d %d", &(DS->Riv[i].shape), &(DS->Riv[i].material));          /*    Read Shape and Material Type          */
        fscanf(riv_file, "%d %d", &(DS->Riv[i].IC), &(DS->Riv[i].BC));                   /*    Read Initial and Boundary Condition   */
        fscanf(riv_file, "%d", &(DS->Riv[i].reservoir));                                 /*    Read Reservoir Type                   */
    }

    fscanf(riv_file, "%s %d", tempchar, &DS->NumRivShape);                               /*    Read Number of Shape Types            */
    DS->Riv_Shape = (river_shape *)malloc(DS->NumRivShape*sizeof(river_shape));          /*    Allocate Memory for Shape Types       */

    for (i=0; i<DS->NumRivShape; i++)
    {
        fscanf(riv_file, "%d %lf", &DS->Riv_Shape[i].index, &DS->Riv_Shape[i].width);    /*    Read Shape # and Width                */
        fscanf(riv_file, "%lf %lf", &DS->Riv_Shape[i].depth, &DS->Riv_Shape[i].bed);     /*    Read Depth and Bed Elevation          */
        fscanf(riv_file, "%d %lf",&DS->Riv_Shape[i].interpOrd,&DS->Riv_Shape[i].coeff);  /*    Read Interpolation Order & Coeff      */

        /*    CALIBRATION    */
        DS->Riv_Shape[i].depth = rivDepth_CALIB*DS->Riv_Shape[i].depth;
        DS->Riv_Shape[i].coeff = rivCoeff_CALIB*DS->Riv_Shape[i].coeff;
    }

    fscanf(riv_file, "%s %d", tempchar, &DS->NumRivMaterial);                            /* Read Number of Material Types            */
    DS->Riv_Mat = (river_material *)malloc(DS->NumRivMaterial*sizeof(river_material));   /* Allocate Memory for Material Types       */

    for (i=0; i<DS->NumRivMaterial; i++)
    {
        fscanf(riv_file, "%d %lf %lf %lf", &DS->Riv_Mat[i].index, &DS->Riv_Mat[i].Rough, &DS->Riv_Mat[i].Cwr, &DS->Riv_Mat[i].Sf);
                                                              /* Read Roughness, Coeff. of Discharge and Manning's Roughness Coeff   */

        /*    CALIBRATION    */
        DS->Riv_Mat[i].Rough=roughRiv_CALIB*DS->Riv_Mat[i].Rough;
    }

    fscanf(riv_file, "%s %d", tempchar, &DS->NumRivIC);                                  /* Read Number of Initial Condition Types   */
    DS->Riv_IC = (river_IC *)malloc(DS->NumRivIC*sizeof(river_IC));                      /* Allocate Memory for Initial Cond. Types  */

    for (i=0; i<DS->NumRivIC; i++)
    {
        fscanf(riv_file, "%d %lf", &DS->Riv_IC[i].index, &DS->Riv_IC[i].value);          /* Read Initial Condition Value             */
    }

    fscanf(riv_file, "%s %d", tempchar, &DS->NumRivBC);                                  /* Read Number of Boundary Conditions       */
    DS->TSD_Riv = (TSD *)malloc(DS->NumRivBC*sizeof(TSD));                               /* Allocated memory for Boundary Conditions */

    for(i=0; i<DS->NumRivBC; i++)
    {
        fscanf(riv_file, "%s %d %d", DS->TSD_Riv[i].name, &DS->TSD_Riv[i].index, &DS->TSD_Riv[i].length);
                                                                                         /* Read Boundary Cond. TimeSeries Lenght    */

        DS->TSD_Riv[i].TS = (realtype **)malloc((DS->TSD_Riv[i].length)*sizeof(realtype));    /* Allocate Memory for TS              */
        for(j=0; j<DS->TSD_Riv[i].length; j++)
        {
            DS->TSD_Riv[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));                    /* Allocate Memory for TS Cont...      */
        }

        for(j=0; j<DS->TSD_Riv[i].length; j++)
        {
            fscanf(riv_file, "%lf %lf", &DS->TSD_Riv[i].TS[j][0], &DS->TSD_Riv[i].TS[j][1]);  /* Read Time & BC Value                */
        }
    }

    /* read in reservoir information    */
    fscanf(riv_file, "%s %d", tempchar, &DS->NumRes);
    if(DS->NumRes > 0)
    {
        /* YET TO BE IMPLIMENTED */
    }
    fclose(riv_file);
    printf("done.\n");
    /* Finish reading .riv File */

    /****************************************/
    /*========== open *.forc file ==========*/
    /****************************************/
    printf("\n  5) reading %s.forc ... ", filename);
    fn[5] = (char *)malloc((strlen(filename)+5)*sizeof(char));
    strcpy(fn[5], filename);
    forc_file = fopen(strcat(fn[5], ".forc"), "r");

    if(forc_file == NULL)
    {
        printf("\n  Fatal Error: %s.forc is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading .forc File */
    fscanf(forc_file, "%d %d", &DS->NumPrep, &DS->NumTemp);                /*    Read Number of Precipitaion & Temperature TS     */
    fscanf(forc_file, "%d %d", &DS->NumHumidity, &DS->NumWindVel);         /*     Read Number of Rel. Humidity & Wind Vel. TS     */
    fscanf(forc_file, "%d %d", &DS->NumRn, &DS->NumG);                     /*    Read Number of Solar Radiation & Dummy           */
    fscanf(forc_file, "%d %d", &DS->NumP, &DS->NumLC);                     /*    Read Number of Pressure and LAI TS               */
    fscanf(forc_file, "%d", &DS->NumMeltF);                                /*     Read Number of Melting Factor TS                */
    fscanf(forc_file, "%d", &DS->NumSource);                               /*    Read Number of Sources/Sinks TS                  */

    DS->TSD_Prep = (TSD *)malloc(DS->NumPrep*sizeof(TSD));                 /*    Allocate Memory for Precipitation Time Series    */
    DS->TSD_Temp = (TSD *)malloc(DS->NumTemp*sizeof(TSD));                 /*    Allocate Memory for Temperature Time Series      */
    DS->TSD_Humidity = (TSD *)malloc(DS->NumHumidity*sizeof(TSD));         /*    Allocate Memory for Rel. Humidity Time Series    */
    DS->TSD_WindVel = (TSD *)malloc(DS->NumWindVel*sizeof(TSD));           /*    Allocate Memory for Wind Velocity Time Series    */
    DS->TSD_Rn = (TSD *)malloc(DS->NumRn*sizeof(TSD));                     /*    Allocate Memory for Solar Radiation Time Series  */
    DS->TSD_G = (TSD *)malloc(DS->NumG*sizeof(TSD));                       /*    Allocate Memory for DUMMY Time Series            */
    DS->TSD_Pressure = (TSD *)malloc(DS->NumP*sizeof(TSD));                /*    Allocate Memory for Pressure Time Series         */
    DS->TSD_LAI = (TSD *)malloc(DS->NumLC*sizeof(TSD));                    /*    Allocate Memory for Leaf Area Index Time Series  */
    DS->TSD_DH = (TSD *)malloc(DS->NumLC*sizeof(TSD));                     /*    Allocate Memory for Stomal Resist Time Series    */
    DS->TSD_MeltF = (TSD *)malloc(DS->NumMeltF*sizeof(TSD));               /*    Allocate Memory for Melting Factor Time Series   */
    DS->TSD_Source = (TSD *)malloc(DS->NumSource*sizeof(TSD));             /*    Allocate Memory for Sources/Sinks Time Series    */

    DS->SIFactor = (realtype *)malloc(DS->NumLC*sizeof(realtype));         /*    Allocate Memory for SI Factor Time Series        */
    DS->WindH = (realtype *)malloc(DS->NumWindVel*sizeof(realtype));       /*    Allocate Memory for WindVel. Height Time Series  */

    /*    Read All the Precipitation Time Series    */
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
        DS->TSD_Prep[i].iCounter=0;
    }

    /*    Read All the Temperature Time Series    */
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
        DS->TSD_Temp[i].iCounter=0;
    }

    /*    Read All the Rel. Humidity Time Series    */
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
            DS->TSD_Humidity[i].TS[j][1] = (DS->TSD_Humidity[i].TS[j][1] > 1.0?1.0:DS->TSD_Humidity[i].TS[j][1]);
        }
        DS->TSD_Humidity[i].iCounter=0;
    }

    /*    Read All the Wind Velocity Time Series    */
    for(i=0; i<DS->NumWindVel; i++)
    {
        fscanf(forc_file, "%s %d %d %lf", DS->TSD_WindVel[i].name, &DS->TSD_WindVel[i].index, &DS->TSD_WindVel[i].length, &DS->WindH[i]);

        DS->TSD_WindVel[i].TS = (realtype **)malloc((DS->TSD_WindVel[i].length)*sizeof(realtype));

        for(j=0; j<DS->TSD_WindVel[i].length; j++)
        {
            DS->TSD_WindVel[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
        }

        for(j=0; j<DS->TSD_WindVel[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_WindVel[i].TS[j][0], &DS->TSD_WindVel[i].TS[j][1]);
        }
        DS->TSD_WindVel[i].iCounter=0;
    }

    /*    Read All the Solar Radiation Time Series    */
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
        DS->TSD_Rn[i].iCounter=0;
    }

    /*    Read All the DUMMY Time Series    */
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
        DS->TSD_G[i].iCounter=0;
    }

    /*    Read All the Vapor Pressure Time Series    */
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
        DS->TSD_Pressure[i].iCounter=0;
    }

    /*    Read All the LAI Time Series    */
    for(i=0; i<DS->NumLC; i++)
    {
        fscanf(forc_file, "%s %d %d %lf", DS->TSD_LAI[i].name, &DS->TSD_LAI[i].index, &DS->TSD_LAI[i].length, &DS->SIFactor[i]);
        DS->TSD_LAI[i].TS = (realtype **)malloc((DS->TSD_LAI[i].length)*sizeof(realtype));

        for(j=0; j<DS->TSD_LAI[i].length; j++)
        {
            DS->TSD_LAI[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
        }

        for(j=0; j<DS->TSD_LAI[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_LAI[i].TS[j][0], &DS->TSD_LAI[i].TS[j][1]);

            DS->TSD_LAI[i].TS[j][1]=lai_CALIB*DS->TSD_LAI[i].TS[j][1];
        }
        DS->TSD_LAI[i].iCounter=0;
    }

    /*    Read All the Displacement Height Time Series    */
    for(i=0; i<DS->NumLC; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_DH[i].name, &DS->TSD_DH[i].index, &DS->TSD_DH[i].length);

        DS->TSD_DH[i].TS = (realtype **)malloc((DS->TSD_DH[i].length)*sizeof(realtype));

        for(j=0; j<DS->TSD_DH[i].length; j++)
        {
            DS->TSD_DH[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
        }
        for(j=0; j<DS->TSD_DH[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_DH[i].TS[j][0], &DS->TSD_DH[i].TS[j][1]);
        }
        DS->TSD_DH[i].iCounter=0;
    }

    /*    Read All the Melting Factor Time Series    */
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
        DS->TSD_MeltF[i].iCounter=0;
    }

    /*    Read All the Sources/Sinks Time Series    */
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

        DS->TSD_Source[i].iCounter=0;
    }

    fclose(forc_file);
    printf("done.\n");
    /* Finish reading .forc File */

    /***************************************/
    /*========== open *.ibc file ==========*/
    /***************************************/
    printf("\n  6) reading %s.ibc  ... ", filename);
    fn[6] = (char *)malloc((strlen(filename)+4)*sizeof(char));
    strcpy(fn[6], filename);
    ibc_file =  fopen(strcat(fn[6], ".ibc"), "r");

    if(ibc_file == NULL)
    {
        printf("\n  Fatal Error: %s.ibc is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading .ibc File */
    fscanf(ibc_file, "%d %d", &DS->Num1BC, &DS->Num2BC);

    if(DS->Num1BC+DS->Num2BC > 0)
    {
        DS->TSD_EleBC = (TSD *)malloc((DS->Num1BC+DS->Num2BC)*sizeof(TSD));
    }

    if(DS->Num1BC>0)
    {    /* For elements with Dirichilet Boundary Conditions */
        for(i=0; i<DS->Num1BC; i++)
        {
            fscanf(ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index,
                                   &DS->TSD_EleBC[i].length);

            DS->TSD_EleBC[i].TS = (realtype **)malloc((DS->TSD_EleBC[i].length)*sizeof(realtype));

            for(j=0; j<DS->TSD_EleBC[i].length; j++)
            {
                DS->TSD_EleBC[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
            }

            for(j=0; j<DS->TSD_EleBC[i].length; j++)
            {
                fscanf(ibc_file, "%lf %lf", &DS->TSD_EleBC[i].TS[j][0],
                                   &DS->TSD_EleBC[i].TS[j][1]);
            }
            DS->TSD_EleBC[i].iCounter=0;
        }
    }

    if(DS->Num2BC>0)
    {    /* For elements with Nueman (non-natural) Boundary Conditions */
        /* This part of code has not be tested ! */
        for(i=DS->Num1BC; i<DS->Num1BC+DS->Num2BC; i++)
        {
            fscanf(ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index,
                                   &DS->TSD_EleBC[i].length);

            DS->TSD_EleBC[i].TS = (realtype **)malloc((DS->TSD_EleBC[i].length)*sizeof(realtype));
            for(j=0; j<DS->TSD_EleBC[i].length; j++)
            {
                DS->TSD_EleBC[i].TS[j] = (realtype *)malloc(2*sizeof(realtype));
            }
            for(j=0; j<DS->TSD_EleBC[i].length; j++)
            {
                fscanf(forc_file, "%lf %lf", &DS->TSD_EleBC[i].TS[j][0], &DS->TSD_EleBC[i].TS[j][1]);
            }
            DS->TSD_EleBC[i].iCounter=0;
        }
    }
    fscanf(ibc_file, "%d", &DS->NumEleIC);
    /*DS->Ele_IC = (element_IC *)malloc(DS->NumEleIC*sizeof(element_IC));
    for(i=0; i<DS->NumEleIC; i++)
    {
        fscanf(ibc_file, "%d", &DS->Ele_IC[i].index);
        fscanf(ibc_file, "%lf", &DS->Ele_IC[i].interception);
        fscanf(ibc_file, "%lf", &DS->Ele_IC[i].snow);
        fscanf(ibc_file, "%lf", &DS->Ele_IC[i].surf);
        fscanf(ibc_file, "%lf", &DS->Ele_IC[i].unsat);
        fscanf(ibc_file, "%lf", &DS->Ele_IC[i].sat);
    }*/

    fclose(ibc_file);
    printf("done.\n");
    /* Finish reading .ibc File */

    /****************************************/
    /*========== open *.para file ==========*/
    /****************************************/
    printf("\n  7) reading %s.para ... ", filename);
    fn[7] = (char *)malloc((strlen(filename)+5)*sizeof(char));
    strcpy(fn[7], filename);
    para_file = fopen(strcat(fn[7], ".para"), "r");

    if(para_file == NULL)
    {
        printf("\n  Fatal Error: %s.para is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading .para File */
    fscanf(para_file, "%d %d", &CS->Verbose, &CS->Debug);
    fscanf(para_file, "%d", &CS->int_type);
    //fscanf(para_file, "%d %d %d %d", &CS->res_out, &CS->flux_out, &CS->q_out, &CS->etis_out);
    fscanf(para_file, "%d %d %d", &DS->UnsatMode, &DS->SurfMode, &DS->RivMode);
    fscanf(para_file, "%d", &CS->Solver);
    if(CS->Solver == 2) //TODO : Needs Correction !!!
    {
        fscanf(para_file, "%d %d %lf", &CS->GSType, &CS->MaxK, &CS->delt);
    }
    fscanf(para_file, "%lf %lf", &CS->abstol, &CS->reltol);
    fscanf(para_file, "%lf %lf %lf", &CS->InitStep, &CS->MaxStep, &CS->ETStep);
    fscanf(para_file, "%lf %lf %d", &CS->StartTime, &CS->EndTime, &CS->outtype);
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
    /* Finish reading .para File */

}

