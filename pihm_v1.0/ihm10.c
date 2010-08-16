/*******************************************************************************
 * File        : ihm10.c                                                       *
 * Function    : This is main entrance of IHM 1.0                              *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               * 
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * IHM1.0 is an implementation of Yizhong Qu's PhD Thesis. It is an integrated *
 * finite volume hydrologic model. It simulates channel routing, overland flow *
 * and groundwater flow in full coupled scheme. It uses semi-discrete approach *
 * to discretize PDE into ODE, and solved it with SUNDIAL (LLNL) package.      *
 *                                                                             *
 * Reference                                                                   *
 *                                                                             *
 * 1.Yizhong Qu, An Integrated Hydrological Model Using Semi-Discrete Finite   *
 *               Volume Formulation, PhD Thesis, the Pennsylvanian State       *
 *               University, 2004                                              *
 * 2.Yizhong Qu, Christopher Duffy, Toward An Integrate Hydrological Model For *
 *               River Basins: A Semi-Discrete, Finite Volume Formulation      *
 *                                                                             *
 * This code is free for users with research purpose only, if appropriate      *
 * citation is refered. However, there is no warranty in any format for this   *
 * product.                                                                    *
 *                                                                             *
 * For questions or comments, please contact the authors of the reference.     *
 * One who want to use it for other consideration may also contact Dr.Duffy    *
 * at cxd11@psu.edu.                                                           *
 *                                                                             *      
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* SUNDIAL Header Files */
#include "sundialstypes.h"   /* realtype, integertype, booleantype defination */
#include "cvode.h"           /* CVODE header file                             */
#include "iterativ.h"        /* for types of preconditioning and Gram-Schmidt */
#include "cvspgmr.h"         /* CVSPGMR linear header file                    */
#include "smalldense.h"      /* use generic DENSE linear solver for "small"   */
#include "nvector_serial.h"  /* contains the definition of type N_Vector      */
#include "sundialsmath.h"    /* contains UnitRoundoff, RSqrt, SQR functions   */
#include "cvdense.h"         /* CVDENSE header file                           */
#include "dense.h"           /* generic dense solver header file              */

/* IHM Header Files */
#include "ihm.h"             /* Definations for all date Structure in IHM     */

/* Function to calculate right hand side of ODE systems */
void initialize(char *, Model_Data, Control_Data *, N_Vector);
void calET(realtype, realtype, N_Vector, Model_Data);
void calIS(realtype, realtype, Model_Data);
void f(integertype, realtype, N_Vector, N_Vector, void *);
void read_alloc(char *, Model_Data, Control_Data *);
  
/* Main Function */
int main(int argc, char *argv[])
{  
  char *filename = "shalehills";        /* Input file name prefix    */
  char *StateFile;                /* Output file name string   */
  char *FluxFile;
  char *ETISFile;  
  char *QFile;  
  
  Model_Data mData;               /* Model Data                */
  Control_Data cData;             /* Solver Control Data       */
  
  N_Vector CV_Y;                  /* State Variables Vector    */
  M_Env machEnv;                  /* Machine Environments      */
  
  realtype ropt[OPT_SIZE];        /* Optional real type and integer type  */
  long int iopt[OPT_SIZE];        /* vecter for Message Passing to solver */
  
  void *cvode_mem;                /* Model Data Pointer        */
  int flag;                       /* flag to test return value */
  
  FILE *res_state_file;           /* Output file for States    */
  FILE *res_flux_file;            /* Output file for Flux      */
  FILE *res_etis_file;            /* Output file for ET and IS */
  FILE *res_q_file;
  
  int N;                          /* Problem size              */
  int i;                          /* loop index                */
  realtype t;                     /* simulation time           */
  realtype NextPtr, StepSize;     /* stress period & step size */
  
  clock_t start, end_r, end_s;    /* system clock at points    */
  realtype cputime_r, cputime_s;  /* for duration in realtype  */
  
  /* allocate memory for model data structure */
  mData = (Model_Data)malloc(sizeof *mData);
  start = clock();
  
  /* get user specified file name in command line */
  if(argc >= 2)
  {
    filename = (char *)malloc(strlen(argv[1])*sizeof(char));
    strcpy(filename, argv[1]);
  }  
  
  printf("\nBelt up!  PIHM 1.0 is starting ... \n");
  
  /* read in 7 input files with "filename" as prefix */
  read_alloc(filename, mData, &cData); 
    
  /* problem size */
  N = 3*mData->NumEle + mData->NumRiv;
  
  /* initial machine environment variable */
  machEnv = M_EnvInit_Serial(N);
  
  /* initial state variable depending on machine*/
  CV_Y = N_VNew(N, machEnv);
  
  /* initialize mode data structure */
  initialize(filename, mData, &cData, CV_Y);
  
  if(cData.Debug == 1) {PrintModelData(mData);}  
  
  end_r = clock();
  cputime_r = (end_r - start)/(realtype)CLOCKS_PER_SEC;
  
  printf("\nSolving ODE system ... \n");
  
  /* initial control parameter for CVODE solver. Otherwise the default value by C could cause problems. */
  for(i=0; i<OPT_SIZE; i++)
  {
    ropt[i] = 0.0;
    iopt[i] = 0;
  }
  
  /* set user specified control parameter */
  ropt[H0] = cData.InitStep;  
  ropt[HMAX] = cData.MaxStep; 
  
  /* allocate memory for solver */
  cvode_mem = CVodeMalloc(N, f, cData.StartTime, CV_Y, BDF, NEWTON, SS, &cData.reltol, 
                          &cData.abstol, mData, NULL, TRUE, iopt, ropt, machEnv);
  if(cvode_mem == NULL) {printf("CVodeMalloc failed. \n"); return(1);}
  
  if(cData.Solver == 1)
  {
    /* using dense direct solver */
    flag = CVDense(cvode_mem, NULL, NULL);
    if(flag != SUCCESS) {printf("CVDense failed. \n"); return(1);}
  } 
  else if(cData.Solver == 2)
  {
    /* using iterative solver */
    flag = CVSpgmr(cvode_mem, NONE, cData.GSType, cData.MaxK, cData.delt, NULL, NULL,
                       mData, NULL, NULL);
    if (flag != SUCCESS) {printf("CVSpgmr failed."); return(1);}
  } 
  
  /*allocate and copy to get output file name */
  StateFile = (char *)malloc((strlen(filename)+4)*sizeof(char));
  strcpy(StateFile, filename);
  FluxFile  = (char *)malloc((strlen(filename)+5)*sizeof(char));
  strcpy(FluxFile, filename);
  ETISFile  = (char *)malloc((strlen(filename)+5)*sizeof(char));
  strcpy(ETISFile, filename);
  QFile = (char *)malloc((strlen(filename)+2)*sizeof(char));
  strcpy(QFile, filename);
  
  /* open output file */
  if (cData.res_out == 1) {res_state_file = fopen(strcat(StateFile, ".res"), "w");}
  if (cData.flux_out == 1) {res_flux_file = fopen(strcat(FluxFile, ".flux"), "w");}
  if (cData.etis_out == 1) {res_etis_file = fopen(strcat(ETISFile, ".etis"), "w");}
  if (cData.q_out == 1) {res_q_file = fopen(strcat(QFile, ".q"), "w");}
  
  /* print header of output file */
  if (cData.res_out == 1) {FPrintYheader(res_state_file, mData);}
  if (cData.etis_out == 1) {FPrintETISheader(res_etis_file, mData);}
  if (cData.q_out == 1) {FPrintETISheader(res_q_file, mData);}
  printf("\n");
  
  /* set start time */
  t = cData.StartTime;
  
  /* start solver in loops */
  for(i=0; i<cData.NumSteps; i++)
  {
    /* prompt information in non-verbose mode */
    if (cData.Verbose != 1)
    {
      printf("  Running: %-4.1f%% ... ", (100*(i+1)/((realtype) cData.NumSteps))); 
      fflush(stdout);
    }
    
    /* inner loops to next output points with ET step size control */
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
      
      /* calculate Interception Storage */
      calIS(t, StepSize, mData);
      /* solving ODE system */
      flag = CVode(cvode_mem, NextPtr, CV_Y, &t, NORMAL);  
      /* calculate ET and adjust states variables*/
      calET(t, StepSize, CV_Y, mData);
    }  
    
    if(cData.Verbose == 1) {PrintVerbose(i, t, iopt, ropt);}

    /* uncomment it if user need it verbose mode */  
    /* if(cData.Verbose == 1) {PrintY(mData, CV_Y, t);} */
    
    /* print out results to files at every output time */
    if (cData.res_out == 1) {FPrintY(mData, CV_Y, t, res_state_file);}
    if (cData.flux_out == 1) {FPrintFlux(mData, t, res_flux_file);}
    if (cData.etis_out == 1) {FPrintETIS(mData, t, res_etis_file);}
    if (cData.q_out == 1) {FPrintQ(mData, t, res_q_file);}
    
    if (cData.Verbose != 1) {printf("\r");}  
    
    if(flag != SUCCESS) {printf("CVode failed, flag = %d. \n", flag); return(flag);} 
    
    /* clear buffer */
    fflush(stdout);  
  }
  
  /* free memory */
  /*
  N_VFree(CV_Y);
  CVodeFree(cvode_mem);
  M_EnvFree_Serial(machEnv); 
  */
  
  /* capture time */
  end_s = clock();
  cputime_s = (end_s - end_r)/(realtype)CLOCKS_PER_SEC;
  
  /* print out simulation statistics */
  PrintFarewell(cData, iopt, ropt, cputime_r, cputime_s);
  if (cData.res_out == 1) {FPrintFarewell(cData, res_state_file, iopt, ropt, cputime_r, cputime_s);}
  
  /* close output files */
  if (cData.res_out == 1)  {fclose(res_state_file);}
  if (cData.flux_out == 1) {fclose(res_flux_file);}
  if (cData.etis_out == 1) {fclose(res_etis_file);}
  if (cData.q_out == 1)    {fclose(res_q_file);}
  
  free(mData);
  
  return 0;
}

