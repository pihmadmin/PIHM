/*******************************************************************************
 * File        : print_function.c                                              *
 * Function    : print out model result to stdout or output files              *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               * 
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * In this file, several print function defined serves to print our model      *
 * result to stdout(monitor) or output files with the same prefix as in input  *
 * files. Users can modify all the following print function as need accordingly*
 * There are 3 files so far: *.res, *.flux, and *.etis. In res file, all state *
 * variables, including Overland depth, Unsaturated and Saturated Storage for  *
 * every element, and River stage for each river segment are saved. All fluxes *
 * are saved in flux file. etis file has all envapotranspiration and inter-    *
 * ception storages in it.                                                     *
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

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "ihm.h"  
#include "cvode.h" 
#include "cvdense.h"
   

void PrintModelData(Model_Data DS)
{
  PrintEle(DS);
  PrintEleAtt(DS);
  PrintNode(DS);
  PrintRiv(DS);
  PrintSoil(DS);
  PrintForcing(DS);
  PrintTS(DS->TSD_Inc, DS->NumInc);
  fflush(stdout);
}

void PrintEle(Model_Data DS)
{
  int i;
  
  printf("\nElements Information: \n");
  printf("\n   Index     N-1     N-2     N-3    Nr-1    Nr-2    Nr-3");
  printf("      Edge-1      Edge-2      Edge-3           X           Y       Z_MIN       Z_MAX        AREA\n");
  printf("\n");
  
  for (i=0; i<DS->NumEle; i++)
  {
    printf("%8d", DS->Ele[i].index);
    printf("%8d%8d%8d", DS->Ele[i].node[0], DS->Ele[i].node[1], DS->Ele[i].node[2]);
    printf("%8d%8d%8d", DS->Ele[i].nabr[0], DS->Ele[i].nabr[1], DS->Ele[i].nabr[2]);

    printf("%12.5f%12.5f%12.5f", DS->Ele[i].edge[0], DS->Ele[i].edge[1], 
                              DS->Ele[i].edge[2]);
    printf("%12.5f%12.5f%12.5f%12.5f", DS->Ele[i].x, DS->Ele[i].y, DS->Ele[i].zmin, 
                                   DS->Ele[i].zmax);
    printf("%12.5f", DS->Ele[i].area);
    printf("\n");
  }
  printf("\n");
}

void PrintEleAtt(Model_Data DS)
{
  int i;
  
  printf("\nElements Attibute Information: \n");
  printf("\nIndex  Soil   LAI    IC    BC  Prep  Temp Humid  WVel    Rn     G     P\n\n");
  
  for (i=0; i<DS->NumEle; i++)
  {  
    printf("%6d%6d%6d", i+1, DS->Ele[i].soil, DS->Ele[i].LAI);
    printf("%6d%6d", DS->Ele[i].IC, DS->Ele[i].BC);
    printf("%6d%6d", DS->Ele[i].prep, DS->Ele[i].temp);
    printf("%6d%6d", DS->Ele[i].humidity, DS->Ele[i].WindVel);
    printf("%6d%6d", DS->Ele[i].Rn, DS->Ele[i].G);
    printf("%6d\n", DS->Ele[i].pressure);
  }  
}


void PrintNode(Model_Data DS)
{
  int i;
  
  printf("\nNode Information: \n");
  printf("\n  Index       X       Y   Z_min   Z_Max\n");
  printf("\n");
  
  for (i=0; i<DS->NumNode; i++)
  {
    printf("%7d", DS->Node[i].index);
    printf("%8.3lf%8.3lf", DS->Node[i].x, DS->Node[i].y);
    printf("%8.3lf%8.3lf", DS->Node[i].zmin, DS->Node[i].zmax);
    printf("\n");
  }
  printf("\n");
}


void PrintSoil(Model_Data DS)
{
  int i;
  
  printf("\nSoil Information: \n");
  printf("\n  Index    Ksat   SitaS   SitaR   Alpha    Beta");
  printf("      Sf Inc_type\n\n");
  
  for(i=0; i<DS->NumSoil; i++)
  {
    printf("%7d", DS->Soil[i].index);
    printf("%8.4f%8.4f%8.4f", DS->Soil[i].Ksat, DS->Soil[i].SitaS, DS->Soil[i].SitaR);
    printf("%8.4f%8.4f%8.4f", DS->Soil[i].Alpha, DS->Soil[i].Beta, DS->Soil[i].Sf);
    printf("%9d", DS->Soil[i].Inf);
    printf("\n");
  }
}


void PrintTS(TSD *Data, int NumTS)
{
  int i, j;
  
  printf("\nTime Series Data Information: \n\n");
  
  for(i=0; i<NumTS; i++)
  {
    printf("%5s%6d%6d\n\n", Data[i].name, Data[i].index, Data[i].length);
    for(j=0; j<Data[i].length; j++)
    {
      printf("%16.6f %16.6f\n", Data[i].TS[j][0], Data[i].TS[j][1]);
    }
    printf("\n");
  }
}


void PrintRiv(Model_Data DS)
{
  int i;
  
  printf("\nRiver Segments Information: \n");
  printf("\n  Index           X           Y           Z");
  printf("       Depth      Length F_node T_node Down L_ele R_ele Shape Mat  IC  BC RES\n");
  printf("\n");
  
  for(i=0; i<DS->NumRiv; i++)
  {
    printf("%7d%12.4f%12.4f%12.4f%12.4f%12.4f", DS->Riv[i].index, DS->Riv[i].x, 
                                                DS->Riv[i].y, DS->Riv[i].zmin, 
                                                DS->Riv[i].depth, DS->Riv[i].Length);
    printf("%7d%7d%5d", DS->Riv[i].FromNode, DS->Riv[i].ToNode, DS->Riv[i].down);
    printf("%6d%6d%6d", DS->Riv[i].LeftEle, DS->Riv[i].RightEle, DS->Riv[i].shape);
    printf("%4d%4d%4d%4d", DS->Riv[i].material, DS->Riv[i].IC, DS->Riv[i].BC, DS->Riv[i].reservoir);
    printf("\n");
  }
}


void PrintForcing(Model_Data DS)
{
  PrintTS(DS->TSD_Prep, DS->NumPrep);
  PrintTS(DS->TSD_Temp, DS->NumTemp);
  PrintTS(DS->TSD_Humidity, DS->NumHumidity);
  PrintTS(DS->TSD_WindVel, DS->NumWindVel);
  PrintTS(DS->TSD_Rn, DS->NumRn);
  PrintTS(DS->TSD_G, DS->NumG);
  PrintTS(DS->TSD_LAI, DS->NumLAI);
}


void PrintDY(Model_Data DS, N_Vector CV_Y, N_Vector CV_Ydot)
{
  int i, index, Num;
  
  Num = 3*DS->NumEle + DS->NumRiv;
  printf("\nY & DY [1 : %d] = \n\n", Num);
  
  printf("Elements Information\n\n");
  printf(" Index        Surf       Unsat         Sat        Surf       Unsat         Sat\n\n");
  
  for(i=0; i<DS->NumEle; i++)
  {
    index = i-(i/DS->NumEle)*DS->NumEle + 1;
    printf("%6d%12.4f%12.4f%12.4f", index, NV_Ith_S(CV_Y, i), NV_Ith_S(CV_Y, i + DS->NumEle), 
                                           NV_Ith_S(CV_Y, i + 2*DS->NumEle));
    printf("%12.4f%12.4f%12.4f\n", NV_Ith_S(CV_Ydot, i), NV_Ith_S(CV_Ydot, i + DS->NumEle), 
                                   NV_Ith_S(CV_Ydot, i + 2*DS->NumEle));
  }  
  
  printf("\n");
  printf("River Segments Information\n\n");
  printf(" Index           Y          DY\n\n");
  for(i=0; i<DS->NumRiv; i++)
  {
    printf("%6d%12.4f%12.4f\n", i+1, NV_Ith_S(CV_Y, i + 3*DS->NumEle), 
                                     NV_Ith_S(CV_Ydot, i + 3*DS->NumEle));
  }
  
}


void PrintY(Model_Data DS, N_Vector CV_Y, realtype t)
{
  int i, index, Num;
  
  Num = 3*DS->NumEle + DS->NumRiv;
  printf("\nt = %8.2f   Y [1 : %d] = \n\n", t, Num);
  
  printf("Elements Information\n\n");
  printf(" Index        Surf       Unsat         Sat\n\n");
  
  for(i=0; i<DS->NumEle; i++)
  {
    index = i-(i/DS->NumEle)*DS->NumEle + 1;
    printf("%6d%12.4f%12.4f%12.4f\n", index, NV_Ith_S(CV_Y, i), 
                                      NV_Ith_S(CV_Y, i + DS->NumEle), 
                                      NV_Ith_S(CV_Y, i + 2*DS->NumEle));
  }  
  
  printf("\n");
  printf("River Segments Information\n\n");
  printf(" Index           Y\n\n");
  for(i=0; i<DS->NumRiv; i++)
  {
    printf("%6d%12.4f\n", i+1, NV_Ith_S(CV_Y, i + 3*DS->NumEle));
  }
  
}

void FPrintYheader(FILE *res_file, Model_Data mData)
{
  int N;
  
  N = 3*mData->NumEle + mData->NumRiv;
  fprintf(res_file, "\nPIHM 1.0 State Varibles Result: \n");
  fprintf(res_file, "\nNumEle = %8d  NumRiv = %8d\n", mData->NumEle, mData->NumRiv);
  fprintf(res_file, "Problem Size N = %8d\n", N); 
  fprintf(res_file, "\n");
}

void FPrintY(Model_Data DS, N_Vector CV_Y, realtype t, FILE *res_file)
{
  int i, Num;
  
  Num = 3*DS->NumEle + DS->NumRiv;

  fprintf(res_file, "Current time = %10.4f", t);
  
  if (DS->NumEle > 0)
  {
    fprintf(res_file, "\n\n  Overland Flow Depth (1 :%8d):\n", DS->NumEle);
    fprintf(res_file, "\n       ");
    for (i=0; i<10; i++)
    {
      fprintf(res_file, "%16d", i+1);
    }  
    fprintf(res_file, "\n       ");
    
    for(i=0; i<DS->NumEle; i++)
    {
      if (i%10 == 0) {fprintf(res_file, "\n %5d ", i/10);}
      fprintf(res_file, "%16.8f", NV_Ith_S(CV_Y, i));
    }
  
    fprintf(res_file, "\n\n  Unsaturated Soil Moisture Equivelant Depth (1 :%8d):\n", DS->NumEle);
    fprintf(res_file, "\n       ");
    for (i=0; i<10; i++)
    {
      fprintf(res_file, "%16d", i+1);
    }  
    fprintf(res_file, "\n       ");
    
    for(i=0; i<DS->NumEle; i++)
    {
      if (i%10 == 0) {fprintf(res_file, "\n %5d ", i/10);}
      fprintf(res_file, "%16.8f", NV_Ith_S(CV_Y, i+DS->NumEle));
    }
  
    fprintf(res_file, "\n\n  Saturated Groundwater Depth (1 :%8d):\n", DS->NumEle);
    fprintf(res_file, "\n       ");
    for (i=0; i<10; i++)
    {
      fprintf(res_file, "%16d", i+1);
    }  
    fprintf(res_file, "\n       ");
    
    for(i=0; i<DS->NumEle; i++)
    {
      if (i%10 == 0) {fprintf(res_file, "\n %5d ", i/10);}
      fprintf(res_file, "%16.8f", NV_Ith_S(CV_Y, i+2*DS->NumEle));
    }
  }
  
  if (DS->NumRiv > 0)
  {
    fprintf(res_file, "\n\n  Channel Flow Depth (1 :%8d):\n", DS->NumRiv);
    fprintf(res_file, "\n       ");
    for (i=0; i<10; i++)
    {
      fprintf(res_file, "%16d", i+1);
    }  
    fprintf(res_file, "\n       ");
    
    for(i=0; i<DS->NumRiv; i++)
    {
      if (i%10 == 0) {fprintf(res_file, "\n %5d ", i/10);}
      fprintf(res_file, "%16.8f", NV_Ith_S(CV_Y, i+3*DS->NumEle));
    }
  
    /*fprintf(res_file, "%16.8f", DS->Q);*/
  }  
  fprintf(res_file, "\n\n  Discharge at outlet = %16.8f\n", DS->Q);
  fprintf(res_file, "\n");
  fflush(res_file);
}

void FPrintFlux(Model_Data DS, realtype t, FILE *res_file)
{
  int i, j, Num;
  
  fprintf(res_file, "t = %10.4f\n\n", t);
  Num = DS->NumEle;
  
  if (Num > 0) {fprintf(res_file, "  FluxSurf = \n");}
  for(i=0; i<Num; i++)
  {
    fprintf(res_file, "%6d", i+1); 
    for(j=0; j<3; j++)
    {
      fprintf(res_file, "%16.8f", DS->FluxSurf[i][j]);
    }
    fprintf(res_file, "\n");     
  }
  
  if (Num > 0) {fprintf(res_file, "\n  FluxSub = \n");}
  for(i=0; i<Num; i++)
  {
    fprintf(res_file, "%6d", i+1); 
    for(j=0; j<3; j++)
    {
      fprintf(res_file, "%16.8f", DS->FluxSub[i][j]);
    }
    fprintf(res_file, "\n");     
  }
  
  Num = DS->NumRiv;
  
  if (Num > 0) {fprintf(res_file, "\n  FluxRiv = \n");}
  for(i=0; i<Num; i++)
  {
    fprintf(res_file, "%6d", i+1); 
    for(j=0; j<6; j++)
    {
      fprintf(res_file, "%16.8f", DS->FluxRiv[i][j]);
    }
    fprintf(res_file, "\n");     
  }
  fprintf(res_file, "\nQ at outlet = %16.8f\n", DS->Q);
  fprintf(res_file, "\n");
  fflush(res_file);
}

void FPrintETISheader(FILE *res_file, Model_Data DS)
{
  int Num;
  
  Num = DS->NumEle;
  
  fprintf(res_file, "\nIHM 1.0 ET and Interception Result: \n");
  fprintf(res_file, "\nNumEle = %8d  NumRiv = %8d\n", DS->NumEle, DS->NumRiv);
  fprintf(res_file, "Problem Size N = %8d", 3*Num + DS->NumRiv); 
  fprintf(res_file, "\n");
}

void FPrintETIS(Model_Data DS, realtype t, FILE *res_file)
{
  int i, j, Num;
  
  fprintf(res_file, "\nCurrent time = %10.4f", t);
  Num = DS->NumEle;

  if (Num > 0)
  {
    fprintf(res_file, "\n\n  Interception Storage (1 :%8d):\n", DS->NumEle);
    fprintf(res_file, "\n       ");
    for (i=0; i<10; i++)
    {
      fprintf(res_file, "%16d", i+1);
    }  
    fprintf(res_file, "\n       ");
  
    for(i=0; i<Num; i++)
    {
      if (i%10 == 0) {fprintf(res_file, "\n %5d ", i/10);}
      fprintf(res_file, "%16.10f", DS->EleIS[i]);
    }

    fprintf(res_file, "\n\n  Envapotranpiration (1 :%8d):\n", DS->NumEle);
  
    fprintf(res_file, "\n       ");
    for (i=0; i<10; i++)
    {
      for (j=0; j<4; j++)
      {
        fprintf(res_file, "%14d_%1d", i+1, j);
      }  
    }  
    fprintf(res_file, "\n       ");
  
    for(i=0; i<Num; i++)
    {
      if (i%10 == 0) {fprintf(res_file, "\n %5d ", i/10);}
      for(j=0; j<4; j++)
      {
        fprintf(res_file, "%16.10f", DS->EleET[i][j]);
      }
    }
  }
  fprintf(res_file, "\n");
  fflush(res_file);
}

void FPrintQ(Model_Data DS, realtype t, FILE *res_file)
{
  fprintf(res_file, "\ntime = %10.4f,  Q = %-16.8f", t, DS->Q);
  fflush(res_file);
}

void PrintVerbose(int i, realtype t, long int iopt[], realtype ropt[])
{
  printf("  Step %5d: ", i+1);  
  printf("t =%10.4lf  t_r =%10.4lf", t, ropt[TCUR]);
  printf("  DeltT =%6.4lf  N_DeltT =%6.4lf", ropt[HU], ropt[HCUR]);
  printf("  Order = %1ld",  iopt[QU]);
  printf("  InterStep =%6ld  RHS Eval =%6ld", iopt[NST], iopt[NFE]);
  printf("\n");
}      

void PrintFarewell(Control_Data cData, long int iopt[], realtype ropt[], realtype cputime_r, realtype cputime_s)
{
  if (cData.Verbose != 1) {printf("\n");}
  printf("\nCongratulations! Simulation Finishes Successfully. \n");
  if(cData.Verbose == 1)
  { 
    PrintFinalStats(iopt, ropt); 
  }
  
  if (cData.Verbose != 1) {printf("\n");}
  if(cData.Solver == 1) {printf("Full Dense Direct Solver.\n");}
  else if(cData.Solver == 2) {printf("Iterative GMRES Solver.\n");}
  printf("Time taken to read in files    = %10.4f seconds. \n", cputime_r);
  printf("Time taken to solve ODE system = %10.4f seconds. \n", cputime_s);  
  printf("\nThanks for running PIHM 1.0. \n");
}

/* function to print out Final Status to stdout */
void PrintFinalStats(long int iopt[], realtype ropt[])
{
  printf("\nFinal Statistics.. \n\n");
  printf("Num Int Steps = %-6ld  Num Linear Sol Setup = %-6ld  RHS Eval   = %-6ld  Num.Jac.Eval   = %-6ld\n",
	 iopt[NST], iopt[NSETUPS], iopt[NFE], iopt[DENSE_NJE]);
  printf("Newton Iter # = %-6ld  Nonlinear Failure #  = %-6ld  Err Failed = %-6ld  Last Step Size = %-6.2lf\n",
	 iopt[NNI], iopt[NCFN], iopt[NETF], ropt[HU]);
  printf("Real Workspace Size: CVODE = %-6ld Solver = %-6ld\n", iopt[LENRW], iopt[DENSE_LRW]); 
  printf("Int Workspace Size:  CVODE = %-6ld Solver = %-6ld\n\n", iopt[LENIW], iopt[DENSE_LIW]);
}

void FPrintFarewell(Control_Data cData, FILE *res_file, long int iopt[], realtype ropt[], realtype cputime_r, realtype cputime_s)
{
  fprintf(res_file, "\nCongratulations! Simulation Finishes Successfully. \n");
  FPrintFinalStats(res_file, iopt, ropt); 
  fprintf(res_file, "The Solver is ");
  
  if(cData.Solver == 1) {fprintf(res_file, "Full Dense Direct Solver.\n");}
  else if(cData.Solver == 2) {fprintf(res_file, "Iterative GMRES Solver.\n");}
  
  fprintf(res_file, "Time taken to read in files    = %10.4f seconds. \n", cputime_r);
  fprintf(res_file, "Time taken to solve ODE system = %10.4f seconds. \n", cputime_s);  
  fprintf(res_file, "\nThanks for running PIHM 1.0. \n");  
}

/* function to print out Final Status to output files */
void FPrintFinalStats(FILE *res_file, long int iopt[], realtype ropt[])
{
  fprintf(res_file, "\nFinal Statistics.. \n\n");
  fprintf(res_file, "Num Int Steps = %-6ld  Num Linear Sol Setup = %-6ld  RHS Eval   = %-6ld  Num.Jac.Eval   = %-6ld\n",
	 iopt[NST], iopt[NSETUPS], iopt[NFE], iopt[DENSE_NJE]);
  fprintf(res_file, "Newton Iter # = %-6ld  Nonlinear Failure #  = %-6ld  Err Failed = %-6ld  Last Step Size = %-6.2lf\n",
	 iopt[NNI], iopt[NCFN], iopt[NETF], ropt[HU]);
  fprintf(res_file, "Real Workspace Size: CVODE = %-6ld Solver = %-6ld\n", iopt[LENRW], iopt[DENSE_LRW]); 
  fprintf(res_file, "Int Workspace Size:  CVODE = %-6ld Solver = %-6ld\n\n", iopt[LENIW], iopt[DENSE_LIW]);
}

