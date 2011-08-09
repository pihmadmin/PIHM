/*******************************************************************************
 * File        : initialize.c                                                  *
 * Function    : initialization of data structue                               *
 * Programmers : Yizhong Qu   @ Pennsylvania State Univeristy                  *
 *               Mukesh Kumar @ Pennsylvania State Univeristy                  *
 *               Gopal Bhatt  @ Pennsylvania State Univeristy                  *
 * Version     : 2.0 (July 10, 2007)                                           *
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

//! @file initialize.c dependent variables of model and control data structure are initialized

/* C Header Files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* SUNDIALS Header Files */
#include "sundials_types.h"
#include "nvector_serial.h"

/* PIHM Header Files */
#include "pihm.h"
#include "calib.h"

/* Calibration Parameters */
realtype satD_CALIB;
realtype br_CALIB;
realtype poros_CALIB;
realtype icsat_CALIB;
realtype rivEle_CALIB;

int lbool;    /**< Optional: To find Sinks    */

/*******************************************************************************
*    Initialize Model & Control Data
********************************************************************************/
void initialize(char *filename, Model_Data DS, Control_Data *CS, N_Vector CV_Y)
//! Function initializes several dependent variables of Model and Control Data Structure
/*! \param filename is Identifier of input files
    \param DS is pointer to model data structure
    \param CS is pointer to control data structure
    \param CV_Y	is state variable vector
*/
{
      int i,j,counterMin, counterMax, MINCONST, domcounter;
      realtype a_x, a_y, b_x, b_y, c_x, c_y, MAXCONST;
      realtype a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax;
      realtype tempvalue;
      FILE *int_file;
      char *fn;

      realtype *zmin_cor;

    /* initializing calibration parameters    */
      satD_CALIB = setsatD_CALIB();
      br_CALIB   = setbr_CALIB();
      poros_CALIB = setporos_CALIB();
      icsat_CALIB = seticsat_CALIB();
      rivEle_CALIB = setrivEle_CALIB();

      zmin_cor=(realtype *)malloc(DS->NumEle*sizeof(realtype));

      printf("\nInitializing data structure ... ");

      for(i=0; i<DS->NumEle; i++)
      {
        a_x = DS->Node[DS->Ele[i].node[0]-1].x;          /* x coordinate of node 1       */
        b_x = DS->Node[DS->Ele[i].node[1]-1].x;          /* x coordinate of node 2       */
        c_x = DS->Node[DS->Ele[i].node[2]-1].x;          /* x coordinate of node 3       */
        a_y = DS->Node[DS->Ele[i].node[0]-1].y;          /* y coordinate of node 1       */
        b_y = DS->Node[DS->Ele[i].node[1]-1].y;          /* y coordinate of node 2       */
        c_y = DS->Node[DS->Ele[i].node[2]-1].y;          /* y coordinate of node 3       */

        a_zmin = DS->Node[DS->Ele[i].node[0]-1].zmin;    /* Bed elevation of node 1      */
        b_zmin = DS->Node[DS->Ele[i].node[1]-1].zmin;    /* Bed elevation of node 2      */
        c_zmin = DS->Node[DS->Ele[i].node[2]-1].zmin;    /* Bed elevation of node 3      */
        a_zmax = DS->Node[DS->Ele[i].node[0]-1].zmax;    /* Surface elevation of node 1  */
        b_zmax = DS->Node[DS->Ele[i].node[1]-1].zmax;    /* Surface elevation of node 2  */
        c_zmax = DS->Node[DS->Ele[i].node[2]-1].zmax;    /* Surface elevation of node 3  */

        /* Finding Lowest and Highest Elevation of an Element and Distance between them    */
        MINCONST =  10000000;
        MAXCONST = -1000000;
        for(j=0; j<3; j++)
        {
            if(DS->Node[DS->Ele[i].node[j]-1].zmin<MINCONST)
            {
                MINCONST=DS->Node[DS->Ele[i].node[j]-1].zmin;
                counterMin=j;
            }
            if(DS->Node[DS->Ele[i].node[j]-1].zmax>MAXCONST)
            {
                MAXCONST=DS->Node[DS->Ele[i].node[j]-1].zmax;
                counterMax=j;
            }
        }
        DS->Ele[i].NodeZmin= DS->Node[DS->Ele[i].node[counterMin]-1].zmin;            /*    Lowest Elevation      */
        DS->Ele[i].NodeZmax= DS->Node[DS->Ele[i].node[counterMax]-1].zmax;            /*    Highest Elevation     */
        DS->Ele[i].NodeDist= sqrt(pow(DS->Node[DS->Ele[i].node[counterMin]-1].x-DS->Node[DS->Ele[i].node[counterMax]-1].x,2)+pow(DS->Node[DS->Ele[i].node[counterMin]-1].y-DS->Node[DS->Ele[i].node[counterMax]-1].y,2));
                                                              /*    Distance between Lowest & Highest Ele Nodes    */
        DS->Ele[i].area = 0.5*((b_x - a_x)*(c_y - a_y) - (b_y - a_y)*(c_x - a_x));    /*    Area of the Element    */
        //printf("\n%lf",DS->Ele[i].area);
        DS->Ele[i].zmax = (a_zmax + b_zmax + c_zmax)/3.0;                             /*    Mean Surface Elevation of an element    */
        DS->Ele[i].zmin = (a_zmin + b_zmin + c_zmin)/3.0;                             /*    Mean Bed Elevation of an element        */

         //DS->Ele[i].zmin =DS->Ele[i].zmax-br_CALIB;


        /* Calculate centroid of triangle */
        DS->Ele[i].x = (a_x + b_x + c_x)/3.0;
        DS->Ele[i].y = (a_y + b_y + c_y)/3.0;


        /*    Calculate Edge Lengths    */
        DS->Ele[i].edge[0] = pow((b_x - c_x), 2) + pow((b_y - c_y), 2);                /*    Length of Edge 1 of an element    */
        DS->Ele[i].edge[1] = pow((c_x - a_x), 2) + pow((c_y - a_y), 2);                /*    Length of Edge 2 of an element    */
        DS->Ele[i].edge[2] = pow((a_x - b_x), 2) + pow((a_y - b_y), 2);                /*    Length of Edge 3 of an element    */


        /* calculate circumcenter of triangle */
          /*DS->Ele[i].x = a_x - ((b_y - a_y)*DS->Ele[i].edge[2] - (c_y - a_y)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
        DS->Ele[i].y = a_y + ((b_x - a_x)*DS->Ele[i].edge[2] - (c_x - a_x)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
        */
        DS->Ele[i].edge[0] = sqrt(DS->Ele[i].edge[0]);
        DS->Ele[i].edge[1] = sqrt(DS->Ele[i].edge[1]);
        DS->Ele[i].edge[2] = sqrt(DS->Ele[i].edge[2]);

        /***************** Temporary ET inputs ***********************************/
        DS->Ele[i].windH = DS->WindH[DS->Ele[i].WindVel - 1]; //10;
        /*************************************************************************/

      }
      /*********************************************/
    /*    Optional :: Routine to Find Sinks         */
    /*********************************************/
    /*for(i=0; i<DS->NumEle; i++){
        lbool=0;
        for(j=0; j<3; j++){
            if(DS->Ele[i].nabr[j]>0){
                if(DS->Ele[i].zmax+1.0>DS->Ele[DS->Ele[i].nabr[j]-1].zmax)
                    lbool=1;
            }
        }
        if(lbool==0){
            printf("%d is sink\n",i+1);
        }
    }
    //getchar();
    */

      /*    Memory allocation for flux terms     */
      DS->FluxSurf = (realtype **)malloc(DS->NumEle*sizeof(realtype));       /* Memory allocation for Surface flux                    */
      DS->FluxSub = (realtype **)malloc(DS->NumEle*sizeof(realtype));        /* Memory allocation for Sursurface flux                 */
      DS->FluxRiv = (realtype **)malloc(DS->NumRiv*sizeof(realtype));        /* Memory allocation for River Flux                      */
      DS->EleET = (realtype **)malloc(DS->NumEle*sizeof(realtype));          /* Memory allocation for Evapotranspiration              */

      for(i=0; i<DS->NumEle; i++)
      {
        DS->FluxSurf[i] = (realtype *)malloc(3*sizeof(realtype));            /* Memory allocation for Surface flux: Cont...           */
        DS->FluxSub[i] = (realtype *)malloc(3*sizeof(realtype));             /* Memory allocation for Subsurface flux: Cont...        */
        DS->EleET[i] = (realtype *)malloc(4*sizeof(realtype));               /* Memory allocation for Evapotranspiration: Cont...     */
      }

      for(i=0; i<DS->NumRiv; i++)
      {
        DS->FluxRiv[i] = (realtype *)malloc(6*sizeof(realtype));             /* Memory allocation for River flux: Cont...             */
      }

      DS->ElePrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));         /* Memory allocation for Precipitation to the Element    */
      DS->EleVic = (realtype *)malloc(DS->NumEle*sizeof(realtype));          /* Memory allocation for Infiltration to the Element     */
      DS->Recharge = (realtype *)malloc(DS->NumEle*sizeof(realtype));        /* Memory allocation for Recharge to GW of a kernel      */
      DS->EleIS = (realtype *)malloc(DS->NumEle*sizeof(realtype));           /* Memory allocation for Interception Storage of Element */
      DS->EleISmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));        /* Memory allocation for Maximum Interception Storage    */
      DS->EleSnow = (realtype *)malloc(DS->NumEle*sizeof(realtype));         /* Memory allocation for Snow accumulation of Element    */
      DS->EleTF = (realtype *)malloc(DS->NumEle*sizeof(realtype));           /* Memory allocation for Throughfall of an Element       */
      DS->Ele2IS = (realtype *)malloc(DS->NumEle*sizeof(realtype));          /* Memory allocation for Interception Storage Rate       */
      DS->EleNetPrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));      /* Memory allocation for Net Precipitation to Element    */

      for(i=0; i<DS->NumEle; i++)
      {
        DS->Ele[i].Ksat = DS->Soil[(DS->Ele[i].soil-1)].Ksat;                /* Saturation Hydraulic Conductivity of an Element       */
        DS->Ele[i].Porosity = DS->Soil[(DS->Ele[i].soil-1)].SitaS -
                          DS->Soil[(DS->Ele[i].soil-1)].SitaR;               /* Porosity of an Element                                */
        DS->Ele[i].Porosity = poros_CALIB*DS->Ele[i].Porosity;               /* CALIBRATION                                           */
        DS->Ele[i].Alpha = DS->Soil[(DS->Ele[i].soil-1)].Alpha;              /* Soil Parameter Alpha of an Element                    */
        DS->Ele[i].Beta = DS->Soil[(DS->Ele[i].soil-1)].Beta;                /* Soil Parameter Beta of an Element                     */
        DS->Ele[i].Sf = DS->Soil[(DS->Ele[i].soil-1)].Sf;                    /* Slope Fricition Factor of an Element                  */
        DS->Ele[i].RzD=DS->Soil[DS->Ele[i].soil-1].RzD>(DS->Ele[i].zmax-DS->Ele[i].zmin)?0.5*(DS->Ele[i].zmax-DS->Ele[i].zmin):DS->Soil[DS->Ele[i].soil-1].RzD;  /* Revisit this modification */
                                                                             /* Root Zone Depth of the Element                        */
        DS->Ele[i].LAImax = DS->LandC[DS->Ele[i].LC-1].LAImax;               /* Maximum Leaf Area Index of the element                */
        DS->Ele[i].Rmin = DS->LandC[DS->Ele[i].LC-1].Rmin;                   /* Minimum Stomatal Resistance of the element            */
        DS->Ele[i].Rs_ref = DS->LandC[DS->Ele[i].LC-1].Rs_ref;               /* Reference ??        */
        DS->Ele[i].Albedo = DS->LandC[DS->Ele[i].LC-1].Albedo;               /* Albedo of an element                                  */
        DS->Ele[i].VegFrac = DS->LandC[DS->Ele[i].LC-1].VegFrac;             /* Vegetation Fraction/Cover of an Element               */
        DS->Ele[i].Rough = DS->LandC[DS->Ele[i].LC-1].Rough;                 /* Manning's Roughness Coefficient of an Element         */
      }

      for(i=0; i<DS->NumRiv; i++)
      {
        DS->Riv[i].x = (DS->Node[DS->Riv[i].FromNode-1].x +
                    DS->Node[DS->Riv[i].ToNode-1].x)/2;                      /* x-centroid of a river segment                         */
        DS->Riv[i].y = (DS->Node[DS->Riv[i].FromNode-1].y +
                    DS->Node[DS->Riv[i].ToNode-1].y)/2;                      /* y-centroid of a river segment                         */
        DS->Riv[i].zmax = (DS->Node[DS->Riv[i].FromNode-1].zmax + DS->Node[DS->Riv[i].ToNode-1].zmax)/2;
                                                                             /* Bank Elevation a river segment                        */
        DS->Riv[i].depth = DS->Riv_Shape[DS->Riv[i].shape-1].depth;          /* Depth of a river segment                              */
        DS->Riv[i].zmin = DS->Riv[i].zmax - DS->Riv[i].depth;                /* Bed Elevation of a river segment                      */

        DS->Riv[i].Length = sqrt(pow(DS->Node[DS->Riv[i].FromNode-1].x -
                                 DS->Node[DS->Riv[i].ToNode-1].x, 2) +
                             pow(DS->Node[DS->Riv[i].FromNode-1].y -
                                 DS->Node[DS->Riv[i].ToNode-1].y, 2));       /* Length of a River Segment                             */
      }
      /************************************************/
    /* Optional : Routine to see River Bed Slope    */
      /************************************************/
      /*for(i=0; i<DS->NumRiv; i++)
      {
        if(DS->Riv[i].down>0)
        {
            if(DS->Riv[i].zmin-DS->Riv[DS->Riv[i].down-1].zmin<=0.0)
            {
                printf("\n%d\t%lf",i,DS->Riv[i].zmin-DS->Riv[DS->Riv[i].down-1].zmin);
                //DS->Ele[DS->Riv[i].LeftEle-1].zmax,DS->Ele[DS->Riv[i].RightEle-1].zmax);
            }
        }
    }
    getchar();
    */

      /* Initialize State Variables */
      if (CS->int_type == 0)      /* relax cases */
      {
        for(i=0; i<DS->NumEle; i++)
        {
              DS->EleIS[i] = 0;                             /* Initialize Interception Storage         */
              NV_Ith_S(CV_Y, i) = 0;                        /* Initialize Surface State                */
              NV_Ith_S(CV_Y, i + DS->NumEle) = (1/DS->Ele[i].Alpha)*(1-exp(-DS->Ele[i].Alpha*(0.1+0.05)));
                                                            /* Initialize Unsaturated Zone State       */
              NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele[i].zmax - DS->Ele[i].zmin - 0.1;
                                                            /* Initialize Saturated Zone State         */
        }

        for(i=0; i<DS->NumRiv; i++)
        {
              NV_Ith_S(CV_Y, i + 3*DS->NumEle) = 0;         /* Initialize River State                  */
        }
      }
      /* type mode */
      else if (CS->int_type == 1)
      {
        if(DS->UnsatMode ==1)
        {
            for(i=0; i<DS->NumEle; i++)
            {
                  DS->EleIS[i] = DS->Ele_IC[i].interception;            /* Initialize Interception Storage        */
                  DS->EleSnow[i]=DS->Ele_IC[i].snow;                    /* Initialize Snow Storage                */
                  NV_Ith_S(CV_Y, i) = DS->Ele_IC[i].surf;               /* Initialize Surface State               */
                  NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[i].sat;   /* Initialize SubSurface State            */
            }

            for(i=0; i<DS->NumRiv; i++)
            {
                  NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Riv_IC[DS->Riv[i].IC-1].value;
                                                                        /* Initialize River State                 */
            }
        }
        if(DS->UnsatMode ==2)
        {
            for(i=0; i<DS->NumEle; i++)
            {
                /*DS->Ele_IC[i].sat= (DS->Ele_IC[i].sat>(DS->Ele[i].zmax - DS->Ele[i].zmin+zmin_cor[i]))?(DS->Ele[i].zmax - DS->Ele[i].zmin)-0.01: (DS->Ele_IC[i].sat-zmin_cor[i]>0?DS->Ele_IC[i].sat-zmin_cor[i]:0);*/
                  //DS->Ele_IC[i].sat=  (DS->Ele[i].zmax - DS->Ele[i].zmin)*icsat_CALIB;
                DS->Ele_IC[i].sat= (DS->Ele[i].zmax-(br_CALIB*(1.0-icsat_CALIB)))>DS->Ele[i].zmin?(DS->Ele[i].zmax-DS->Ele[i].zmin-(br_CALIB*(1.0-icsat_CALIB))):(DS->Ele[i].zmax-DS->Ele[i].zmin)/2.0;
                  DS->EleIS[i] = DS->Ele_IC[i].interception;                 /* Initialize Interception Storage         */
                  NV_Ith_S(CV_Y, i) = DS->Ele_IC[i].surf;                    /* Initialize Surface State                */

                  DS->Ele_IC[i].sat=(DS->Ele[i].zmax - DS->Ele[i].zmin)-(DS->Ele_IC[i].sat+satD_CALIB)>0?((DS->Ele_IC[i].sat+satD_CALIB)<0?0.0:(DS->Ele_IC[i].sat+satD_CALIB)):(DS->Ele[i].zmax - DS->Ele[i].zmin)-0.01;

                  NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[i].unsat+(1-exp(-1.0*DS->Ele[i].Alpha*((DS->Ele[i].zmax - DS->Ele[i].zmin)-DS->Ele_IC[i].sat)))/DS->Ele[i].Alpha; /* delete the later part: Just for Juniata */
                                                                             /* Initialize Unsaturated Zone State    */
                  NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele_IC[i].sat;      /* Initialize Unsaturated Zone State    */

                /* If Sat + UnSat > Bed/Soil Thickness    */
                  if ((NV_Ith_S(CV_Y, i + DS->NumEle) + NV_Ith_S(CV_Y, i + 2*DS->NumEle)) >= (DS->Ele[i].zmax - DS->Ele[i].zmin))
                  {
                    NV_Ith_S(CV_Y, i + DS->NumEle) = ((DS->Ele[i].zmax - DS->Ele[i].zmin) - NV_Ith_S(CV_Y, i + 2*DS->NumEle))*0.9;
                       // printf("\n here %d %e %e %e",i,NV_Ith_S(CV_Y, i + DS->NumEle),NV_Ith_S(CV_Y, i + 2*DS->NumEle),(DS->Ele[i].zmax - DS->Ele[i].zmin));
                    if (NV_Ith_S(CV_Y, i + DS->NumEle) < 0)
                    {
                        NV_Ith_S(CV_Y, i + DS->NumEle) = 0;
                    }
                  }
            }

            for(i=0; i<DS->NumRiv; i++)
            {
                  NV_Ith_S(CV_Y, i + 3*DS->NumEle) = DS->Riv_IC[DS->Riv[i].IC-1].value;    /* Initialize River State      */
                  //NV_Ith_S(CV_Y, DS->Riv[i].LeftEle-1 + 2*DS->NumEle) = (DS->Ele[DS->Riv[i].LeftEle-1].zmax - DS->Ele[DS->Riv[i].LeftEle-1].zmin)*rivEle_CALIB;
                  //NV_Ith_S(CV_Y, DS->Riv[i].RightEle-1 + 2*DS->NumEle) = (DS->Ele[DS->Riv[i].RightEle-1].zmax - DS->Ele[DS->Riv[i].RightEle-1].zmin)*rivEle_CALIB;
                  //NV_Ith_S(CV_Y, DS->Riv[i].LeftEle-1 + DS->NumEle) = (1-exp(-1.0*DS->Ele[i].Alpha*((DS->Ele[i].zmax - DS->Ele[i].zmin)-DS->Ele_IC[i].sat)))/DS->Ele[i].Alpha;
                  //NV_Ith_S(CV_Y, DS->Riv[i].RightEle-1 + DS->NumEle) = (1-exp(-1.0*DS->Ele[i].Alpha*((DS->Ele[i].zmax - DS->Ele[i].zmin)-DS->Ele_IC[i].sat)))/DS->Ele[i].Alpha;
            }
        }
      }
      /* hot start mode */
      else if(CS->int_type == 2)
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
      /**************************************************/
      /* Adding routine to initialize state from a file */
      /* with initCondition/States at any given time    */
      /**************************************************/
      else if(CS->int_type == 3)
      {
        /* Open .init File for read    */
        fn = (char *)malloc((strlen(filename)+5)*sizeof(char));
        strcpy(fn, filename);
        int_file = fopen(strcat(fn, ".init"), "r");

        if(int_file == NULL)    /* Exit if couldn't open the file    */
        {
            printf("\n  Fatal Error: %s.int is in use or does not exist!\n", filename);
            exit(1);
        }
        else
        {
            fscanf(int_file, "%lf", &tempvalue);
            /* Check if the expected Start Time matches with that of the .init file    */
            if(tempvalue != CS->StartTime){
                printf("\n  Fatal Error: Initial time in .init file does not match start time\n");
                exit(1);
            }

            if(DS->UnsatMode == 1)
            {
                for(i=0; i<DS->NumEle; i++)
                {
                    fscanf(int_file, "%lf", &tempvalue);      /* Read Interception Storage State    */
                      DS->EleIS[i] = tempvalue;
                      fscanf(int_file, "%lf", &tempvalue);    /* Read Snow Storage State            */
                      DS->EleSnow[i]=tempvalue;
                      fscanf(int_file, "%lf", &tempvalue);    /* Read Surface Flow State            */
                      NV_Ith_S(CV_Y, i) = tempvalue;
                     /*NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[DS->Ele[i].IC-1].unsat; */
                     fscanf(int_file, "%lf", &tempvalue);     /* Read SubSurface Flow State         */
                      NV_Ith_S(CV_Y, i + DS->NumEle) = tempvalue;
                }

                for(i=0; i<DS->NumRiv; i++)
                {
                    fscanf(int_file, "%lf", &tempvalue);      /* Read River Flow State              */
                    NV_Ith_S(CV_Y, i + 2*DS->NumEle) = tempvalue;
                }
            }
            if(DS->UnsatMode == 2)
            {
                for(i=0; i<DS->NumEle; i++)
                {
                    fscanf(int_file, "%lf", &tempvalue);      /* Read Interception Storage State    */
                    DS->EleIS[i] = tempvalue;
                    fscanf(int_file, "%lf", &tempvalue);      /* Read Snow Storage State            */
                    DS->EleSnow[i]=tempvalue;
                    fscanf(int_file, "%lf", &tempvalue);      /* Read Surface Flow State            */
                    NV_Ith_S(CV_Y, i) = tempvalue;
                    fscanf(int_file, "%lf", &tempvalue);      /* Read Unsaturated Zone State        */
                    NV_Ith_S(CV_Y, i + DS->NumEle) = tempvalue;
                    fscanf(int_file, "%lf", &tempvalue);      /* Read Saturated Zone State          */
                    NV_Ith_S(CV_Y, i + 2*DS->NumEle) = tempvalue;
                }

                for(i=0; i<DS->NumRiv; i++)
                {
                    fscanf(int_file, "%lf", &tempvalue);       /* Read River Flow State             */
                    NV_Ith_S(CV_Y, i + 3*DS->NumEle) = tempvalue;
                }
            }
        }
    }

    printf("done.\n");
}

