/*******************************************************************************/
/*                                                                             */
/*                 8888888b. 8888888 888    888 888b     d888                  */
/*                 888   Y88b  888   888    888 8888b   d8888                  */
/*                 888    888  888   888    888 88888b.d88888                  */
/*                 888   d88P  888   8888888888 888Y88888P888                  */
/*                 8888888P"   888   888    888 888 Y888P 888                  */
/*                 888         888   888    888 888  Y8P  888                  */
/*                 888         888   888    888 888   "   888                  */
/*                 888       8888888 888    888 888       888                  */
/*                                                                             */
/* Version     : v2.0 (July 10, 2007)                                          */
/*                                                                             */
/* File        : pihm.c                                                        */
/* Function    : This is main entrance of PIHM                                 */
/* Programmers : Yizhong Qu   @ Pennsylvania State Univeristy                  */
/*               Mukesh Kumar @ Pennsylvania State Univeristy                  */
/*               Gopal Bhatt  @ Pennsylvania State Univeristy                  */
/*-----------------------------------------------------------------------------*/
/*                                                                             */
/* PIHM is an integrated finite volume hydrologic model. It simulates channel  */
/* routing, overland flow and groundwater flow in full coupled scheme. It uses */
/* semi-discrete approach to discretize PDE into ODE, and solved it with       */
/* SUNDIAL (CVODE 2.2.0) package [http://www.llnl.gov/CASC/sundials/]          */
/*                                                                             */
/* References:                                                                 */
/*                                                                             */
/* 1.Yizhong Qu, An Integrated Hydrological Model Using Semi-Discrete Finite   */
/*               Volume Formulation, PhD Thesis, the Pennsylvanian State       */
/*               University, 2004                                              */
/* 2.Yizhong Qu, Christopher Duffy, Toward An Integrate Hydrological Model For */
/*               River Basins: A Semi-Discrete, Finite Volume Formulation      */
/*                                                                             */
/*                                                                             */
/* This code is free for users with research purpose only, if appropriate      */
/* citation is refered. However, there is no warranty in any format for this   */
/* product.                                                                    */
/*                                                                             */
/* For questions or comments, please contact the authors of the reference.     */
/* One who want to use it for other consideration may also contact Dr. Duffy   */
/* at cxd11@psu.edu.                                                           */
/*                                                                             */
/*******************************************************************************/

/** @file pihm.c
\brief This file is the main entrance of PIHM.

The whole PIHM module consists of 11 files including Makefile. The strategy followed here revolves
around the basic setup for using CVODE Solver. In terms of procedural framework the steps those are
accomplished are as follows:

*1. read_alloc.c reads all the input files

*2. initialize.c initializes the model and control data structure variables.

*3. et_is.c computes interception and evaporation from the canopy

*4. f.c computes the rate of change to state variables to form system of ODEs

*5. CVode() solves the system of ODE to output state variables at specified point of times

*6. state variables are passed on to print.c at the end of each timestep

*This process (steps 3 to 6) of marching in time continues until it reaches the end of simulation time period.
*/

/* C Header Files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



/* SUNDIALS Header Files */
#include "sundials_types.h"             /* contains the definition of the type realtype         */
#include "sundials_dense.h"             /* defines the DenseMat type and  accessor macros       */
#include "sundials_smalldense.h"        /* use generic DENSE linear solver for "small"          */
#include "sundials_math.h"              /* contains UnitRoundoff, RSqrt, SQR functions          */
#include "cvode.h"                      /* header file fpr CVODE                                */
#include "cvode_spgmr.h"                /* the Krylov solver SPGMR in the context of CVODE      */
#include "cvode_dense.h"                /* dense direct linear solver in the context of CVODE   */
#include "nvector_serial.h"             /* defines the serial implementation NVECTOR-SERIAL     */



/* PIHM Header Files */
#include "pihm.h"                       /* Definations for all data Structure in PIHM           */
//#include "et_is.h"

/* Function declarations */
void read_alloc(char *, Model_Data, Control_Data *);     /* read input from files :: read_alloc.c                */
N_Vector N_VNew_Serial(int);                             /**< \brief CVODE::Set vector of initial values         */
void initialize(char *, Model_Data, Control_Data *, N_Vector);
                                                         /* Initialize model & Control Data :: initialize.c      */

void* CVodeCreate(int , int); 					         /**< \brief CVODE::Create the CVODE memory block and specify the Solution Method */
int CVodeSetFdata(void *, void *);                       /**< \brief CVODE::Set Data for right-hand side function              */
int CVodeSetInitStep(void *, realtype);                  /**< \brief CVODE::Set Initial step size                              */
int CVodeSetStabLimDet(void *, booleantype);             /**< \brief CVODE::ON/OFF the BDF stability limit detection algorithm */
int CVodeSetMaxStep(void *, realtype);                   /**< \brief CVODE::Specify the maximum absolute value of the step size*/
int CVodeMalloc(void *, CVRhsFn, realtype, N_Vector, int, realtype, void *); /**< \brief CVODE::provide required problem specifications, allocate internal memory for CVODE, and initialize CVODE*/
int CVSpgmr(void *, int, int);                           /**< \brief CVODE::selects the CVSPGMR linear solver                  */
int CVSpilsSetGSType(void *, int);                       /**< \brief CVODE::specifies Gram-Schmidt orthogonalization to be used*/

void calET_IS(realtype, realtype, Model_Data, N_Vector); /* Calculates ET & IS    :: et_is.c                     */
int CVode(void *, realtype, n_Vector, realtype *, int);  /**< \brief CVODE::Advance solution in time             */
int  f(realtype, N_Vector, N_Vector, void *);            /* RHS of system of ODEs :: f.c                         */

void setTSDiCounter(Model_Data mData, realtype t);       /* set the current position (iCounter) of TSD           */
realtype Interpolation(TSD *Data, realtype t);           /* Data Value at time=t from a TimeSeries               */

void FPrintInit(Model_Data);
void FPrint(Model_Data, N_Vector, realtype);
void FPrintInitFile(Model_Data, Control_Data, N_Vector, int);
void FPrintCloseAll(void);






int satEle, ovrEle;

/* Main Function of PIHM */
int main(int argc, char *argv[])
{
    char *filename;                 /* File name prefix for input/output files    */

    Model_Data mData;               /* Model Data                                 */
    Control_Data cData;             /* Control Data                               */
    N_Vector CV_Y;                  /* State Variables Vector                     */

    void *cvode_mem;                /* pointer to the CVODE memory block          */
    int flag;                       /* return value of cvode function calls       */

    int N;                          /* Problem Size  (Numer of ODEs)              */
    int i,j,k;                        /* loop index variables                     */
    realtype t;                     /* simulation time (real time)                */
    realtype NextPtr, StepSize;     /* stress period & step size                  */

    clock_t start, end_r, end_s;    /* system clock at points                     */ //TODO: get rid of it
    realtype cputime_r, cputime_s;  /* for duration in realtype                   */ //TODO: get rid of it

    /***************************
    Next two lines of variable declarations are for printing flow to estuary/BC */
    realtype loc1_bcEle, loc_Avg_Y_Surf, loc_Avg_Y_Sub, loc_Distance, loc_Dif_Y_Sub, loc_Avg_Ksat, loc_Grad_Y_Sub, Sub_Bdd, loc_Dif_Y_Surf, loc_Grad_Y_Surf, Surf_Bdd;
    int loc_i, loc_j;
    /***************************/

    FILE *base2File, *over2File;

    char tmpFileName[100];
    setFileName(tmpFileName);        /* get File Name specified in calib.c        */
    filename = (char *)malloc(sizeof(char)*strlen(tmpFileName));
    strcpy(filename, tmpFileName);

    printf("\nBelt up!  PIHM 2.0 is starting ... \n");
    start = clock();
    /* allocate memory for model data structure */
    mData = (Model_Data)malloc(sizeof *mData);




    base2File=fopen("sc.base2","w");
    over2File=fopen("sc.over2","w");


    /* read the input files with "filename" as prefix */
    read_alloc(filename, mData, &cData);          /* function definition in read_alloc.c    */

    if(mData->UnsatMode ==1) //take off the option 1 from everywhere
    {
          N = 2*mData->NumEle + mData->NumRiv;    /* Set problem dimension                  */
      }
      if(mData->UnsatMode ==2)
    {
        N = 3*mData->NumEle + mData->NumRiv;      /* Set problem dimension                  */
      }

    CV_Y = N_VNew_Serial(N);                      /* Set Vector of initial values           */


    initialize(filename, mData, &cData, CV_Y);    /* initialize mode data structure         */
                                                  /* function definition in initialize.c    */

    FPrintInit(mData);
    //if(cData.Debug == 1) {PrintModelData(mData);}

    end_r = clock();
    cputime_r = (end_r - start)/(realtype)CLOCKS_PER_SEC;

    printf("\nSolving ODE system ... \n");


    /* Create the CVODE memory block and specify the Solution Method */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if(cvode_mem == NULL) { printf("CVodeCreate failed. \n"); return(1); }


    flag = CVodeSetFdata(cvode_mem, mData);                            /* Set Data for right-hand side function                  */
    flag = CVodeSetInitStep(cvode_mem,cData.InitStep);                 /* Set Initial step size                                  */
    flag = CVodeSetStabLimDet(cvode_mem,TRUE);                         /* ON/OFF the BDF stability limit detection algorithm     */
    flag = CVodeSetMaxStep(cvode_mem,cData.MaxStep);                   /* Specify the maximum absolute value of the step size    */
    flag = CVodeMalloc(cvode_mem, f, cData.StartTime, CV_Y, CV_SS, cData.reltol, &cData.abstol);
                                                                       /* provide required problem specifications,
                                                                         allocate internal memory for CVODE, and initialize CVODE*/
    flag = CVSpgmr(cvode_mem, PREC_NONE, 0);                           /* selects the CVSPGMR linear solver                      */
    flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);                   /* specifies Gram-Schmidt orthogonalization to be used    */


    /*allocate and copy to get output file name */



    t = cData.StartTime;                                               /* set "t" to simulation start time                       */

    /* start CVODE solver in loops for NumStep number of times */
    for(i=0; i<cData.NumSteps; i++)
    {
        /*if (cData.Verbose != 1)
        {
            printf("  Running: %-4.1f%% ... ", (100*(i+1)/((realtype) cData.NumSteps)));
            fflush(stdout);
        }*/

        /* inner loops to next output points with ET/IS step size control */
        while(t < cData.Tout[i+1])
        {
            if (t + cData.ETStep >= cData.Tout[i+1])
            {
                NextPtr = cData.Tout[i+1];
            }
            else
            {
                NextPtr = t + cData.ETStep;
            }
            StepSize = NextPtr - t;

            calET_IS(t, StepSize, mData, CV_Y);                        /* Calculate Evaporation/Interception Rates             */


/******************************************************************************************/
    //tempNetPrep = 0.0;
    //tempPrep    = 0.0;
    //tempET0     = 0.0;
    //tempTF      = 0.0;
    /*
    for(tempI=0; tempI<mData->NumEle; tempI++){
        tempPrep += mData->ElePrep[tempI];
        tempNetPrep += mData->EleNetPrep[tempI];
        tempET0 += mData->EleET[tempI][0]*mData->Ele[tempI].VegFrac;
        tempET1 += mData->EleET[tempI][1];
        tempET2 += mData->EleET[tempI][2];
        tempTF  += mData->EleTF[tempI]*mData->Ele[tempI].VegFrac;
    }
    fprintf(prepFile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, tempPrep/mData->NumEle, tempNetPrep/mData->NumEle, tempET0/mData->NumEle, tempET1/mData->NumEle, tempET2/mData->NumEle,tempTF/mData->NumEle);//fflush(prepFile);
*/
/******************************************************************************************/

            Tsteps=t;
            printf("\n Tsteps = %f ",t);

            flag = CVode(cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);    /* Advance solution in time                                */

            setTSDiCounter(mData, t);
            FPrint(mData, CV_Y, t);


/***************************************************************************************************************/
/*  Specially for Rhode :: Estuary Discharge                                                                   */
/***************************************************************************************************************/
    for(loc_i=0; loc_i<mData->NumEle; loc_i++){
        for(loc_j=0; loc_j<3; loc_j++){
            if(mData->Ele[loc_i].BC > 0){         // Dirichlet BC
                loc1_bcEle = Interpolation(&mData->TSD_EleBC[(mData->Ele[loc_i].BC)-1], t);
                if( loc1_bcEle < mData->Ele[loc_i].zmin ){
                    loc_Avg_Y_Surf = NV_Ith_S(CV_Y, loc_i)/2; //(DummyY[i])/2;
                    loc_Avg_Y_Sub  = NV_Ith_S(CV_Y, loc_i+2*mData->NumEle);//(DummyY[i+2*mData->NumEle])/2;
                }
                else if(loc1_bcEle < mData->Ele[loc_i].zmax){
                    loc_Avg_Y_Surf = NV_Ith_S(CV_Y, loc_i)/2;//(DummyY[i])/2;
                    loc_Avg_Y_Sub  = (loc1_bcEle - mData->Ele[loc_i].zmin + NV_Ith_S(CV_Y, loc_i+2*mData->NumEle))/2.0;
                }
                else{
                    loc_Avg_Y_Surf = (NV_Ith_S(CV_Y, loc_i) + loc1_bcEle - mData->Ele[loc_i].zmax)/2;
                    loc_Avg_Y_Sub  = mData->Ele[loc_i].zmax - mData->Ele[loc_i].zmin;
                }
                loc_Distance = sqrt(pow(mData->Ele[loc_i].edge[0]*mData->Ele[loc_i].edge[1]*mData->Ele[loc_i].edge[2]/(4*mData->Ele[loc_i].area), 2) - pow(mData->Ele[loc_i].edge[loc_j]/2, 2));
                //Sub
                loc_Dif_Y_Sub = NV_Ith_S(CV_Y, loc_i+2*mData->NumEle) + mData->Ele[loc_i].zmin - loc1_bcEle;
                loc_Avg_Ksat = mData->Ele[loc_i].Ksat;
                loc_Grad_Y_Sub = loc_Dif_Y_Sub/loc_Distance;
                Sub_Bdd = loc_Avg_Ksat * loc_Grad_Y_Sub * loc_Avg_Y_Sub * mData->Ele[loc_i].edge[loc_j];

                //Surf
                loc_Dif_Y_Surf = NV_Ith_S(CV_Y, loc_i) + mData->Ele[loc_i].zmax - loc1_bcEle;
                loc_Grad_Y_Surf = loc_Dif_Y_Surf / loc_Distance;
                //Surf_Bdd = 0.1*(Grad_Y_Surf>0?1:-1)* (Avg_Y_Sub * mData->Ele[i].edge[loc_j] ) * (pow(pow(Avg_Y_Surf, 1.0/3.0),2)/(mData->Ele[i].Rough)) * sqrt((Grad_Y_Surf>0?1:-1)*Grad_Y_Surf);

                Surf_Bdd = sqrt(NV_Ith_S(CV_Y, loc_i)/loc_Distance)*pow(NV_Ith_S(CV_Y, loc_i),2.0/3.0)*(NV_Ith_S(CV_Y, loc_i)/2)*mData->Ele[loc_i].edge[loc_j]/mData->Ele[loc_i].Rough;

                fprintf(base2File, "%lf\t", Sub_Bdd);
                fprintf(over2File, "%lf\t", Surf_Bdd);
            }
        }
    }
    fprintf(base2File, "\n");
    fprintf(over2File, "\n");

/***************************************************************************************************************/





            satEle=0;
            ovrEle=0;
            for(k=0; k<mData->NumEle;k++)
            {
                if(NV_Ith_S(CV_Y, k+2*mData->NumEle)/(mData->Ele[k].zmax-mData->Ele[k].zmin)>0.99){
                    satEle++;
                }
                if(NV_Ith_S(CV_Y, k)>1E-4){
                    ovrEle++;
                }
            }


            //fprintf(res_flux_file,"\n"); //fflush(res_flux_file);


/**********************************************MASS BALANCE*****************************************************/

        }


        /* print out results to files at every output time */
           //if (cData.res_out == 1) {FPrintY(mData, CV_Y, t, res_state_file);}
        //if (cData.flux_out == 1) {FPrintFlux(mData, t, res_flux_file);}
        //if (cData.etis_out == 1) {FPrintETIS(mData, t, res_etis_file);}
          /*  if (cData.q_out == 1) {FPrintQ(mData, t, res_q_file);} */

        /*    if (cData.Verbose != 1) {printf("\r");}

        if(flag != SUCCESS) {printf("CVode failed, flag = %d. \n", flag); return(flag);}
          */
        /* clear buffer */
           fflush(stdout);
       }

    FPrintInitFile(mData, cData, CV_Y, i);                    /* Routine for .init File : print.c     */
    FPrintCloseAll();








  /* capture time */
    end_s = clock();
    cputime_s = (end_s - end_r)/(realtype)CLOCKS_PER_SEC;

    /* print out simulation statistics */
    /*PrintFarewell(cData, iopt, ropt, cputime_r, cputime_s);
    if (cData.res_out == 1) {FPrintFarewell(cData, res_state_file, iopt, ropt, cputime_r, cputime_s);}
    */
    /* close output files */
    //if (cData.res_out == 1)  {fclose(res_state_file);}
    //if (cData.flux_out == 1) {fclose(res_flux_file);}
    //if (cData.etis_out == 1) {fclose(res_etis_file);}
    //if (cData.q_out == 1)    {fclose(res_q_file);}

    free(mData);

    return 0;
}








void setTSDiCounter(Model_Data mData, realtype t)
//! Function sets the marker of all the time series according to given time t
/*! \param mData is pointer to model data structure
    \param t is time of simulation
*/
{

    int k;

    for(k=0; k<mData->NumPrep; k++)
    {
        while(mData->TSD_Prep[k].iCounter < mData->TSD_Prep[k].length && t/(24.0*60.0) > mData->TSD_Prep[k].TS[mData->TSD_Prep[k].iCounter+1][0]){
            mData->TSD_Prep[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumTemp; k++){
        while(mData->TSD_Temp[k].iCounter < mData->TSD_Temp[k].length && t/(24.0*60.0) > mData->TSD_Temp[k].TS[mData->TSD_Temp[k].iCounter+1][0]){
            mData->TSD_Temp[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumHumidity; k++){
        while(mData->TSD_Humidity[k].iCounter < mData->TSD_Humidity[k].length && t/(24.0*60.0) > mData->TSD_Humidity[k].TS[mData->TSD_Humidity[k].iCounter+1][0]){
            mData->TSD_Humidity[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumWindVel; k++){
        while(mData->TSD_WindVel[k].iCounter < mData->TSD_WindVel[k].length && t/(24.0*60.0) > mData->TSD_WindVel[k].TS[mData->TSD_WindVel[k].iCounter+1][0]){
            mData->TSD_WindVel[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumRn; k++){
        while(mData->TSD_Rn[k].iCounter < mData->TSD_Rn[k].length && t/(24.0*60.0) > mData->TSD_Rn[k].TS[mData->TSD_Rn[k].iCounter+1][0]){
            mData->TSD_Rn[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumG; k++){
        while(mData->TSD_G[k].iCounter < mData->TSD_G[k].length && t/(24.0*60.0) > mData->TSD_G[k].TS[mData->TSD_G[k].iCounter+1][0]){
            mData->TSD_G[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumP; k++){
        while(mData->TSD_Pressure[k].iCounter < mData->TSD_Pressure[k].length && t/(24.0*60.0) > mData->TSD_Pressure[k].TS[mData->TSD_Pressure[k].iCounter+1][0]){
            mData->TSD_Pressure[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumLC; k++){
        while(mData->TSD_LAI[k].iCounter < mData->TSD_LAI[k].length && t/(24.0*60.0) > mData->TSD_LAI[k].TS[mData->TSD_LAI[k].iCounter+1][0]){
            mData->TSD_LAI[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumLC; k++){
        while(mData->TSD_DH[k].iCounter < mData->TSD_DH[k].length && t/(24.0*60.0) > mData->TSD_DH[k].TS[mData->TSD_DH[k].iCounter+1][0]){
            mData->TSD_DH[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumMeltF; k++){
        while(mData->TSD_MeltF[k].iCounter < mData->TSD_MeltF[k].length && t/(24.0*60.0) > mData->TSD_MeltF[k].TS[mData->TSD_MeltF[k].iCounter+1][0]){
            mData->TSD_MeltF[k].iCounter++;
        }
    }
    for(k=0; k<mData->NumSource; k++){
        while(mData->TSD_Source[k].iCounter < mData->TSD_Source[k].length && t/(24.0*60.0) > mData->TSD_Source[k].TS[mData->TSD_Source[k].iCounter+1][0]){
            mData->TSD_Source[k].iCounter++;
        }
    }

    for(k=0; k<mData->Num1BC+mData->Num2BC; k++){
        while(mData->TSD_EleBC[k].iCounter < mData->TSD_EleBC[k].length && t/(24.0*60.0) > mData->TSD_EleBC[k].TS[mData->TSD_EleBC[k].iCounter+1][0]){
            mData->TSD_EleBC[k].iCounter++;
        }
    }

}
