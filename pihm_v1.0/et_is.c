/*******************************************************************************
 * File        : et_is.c                                                       *
 * Function    : for calculation of evapotranspiration and interception        *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               * 
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * Interception process and evapotranspiration are considered as weakly coupled*
 * processes compared with overland, channel routing and groundwater flow pro- *
 * cesses. Therefore, before each time step, interception is calculated and    *
 * deducted. After each time step, envapotranspiration is calculated and all   *
 * state variables are adjusted accordingly. In this file, ET is calculated    *
 * using Penman-Monteith equation. Most parameters are calculated based on     *
 * Crop Evapotranspiration (FAO No.56).                                        *
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

realtype Interpolation(TSD *Data, realtype t);

void calIS(realtype t, realtype stepsize, void *DS)
{
  int i;
  Model_Data MD;
  
  MD = (Model_Data)DS;
  
  for(i=0; i<MD->NumEle; i++)
  {
    MD->ElePrep[i] = Interpolation(&MD->TSD_Prep[MD->Ele[i].prep-1], t);
    MD->EleISmax[i] = MD->SIFactor[MD->Ele[i].LAI-1]*Interpolation(&MD->TSD_LAI[MD->Ele[i].LAI-1], t);
    
    if(MD->EleIS[i] >= MD->EleISmax[i])
    {
      MD->Ele2IS[i] = 0;
    }
    else if((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + MD->ElePrep[i]*stepsize) >= MD->EleISmax[i]))
    {
      MD->Ele2IS[i] =  (MD->EleISmax[i] - MD->EleIS[i])/stepsize;
      MD->EleIS[i] = MD->EleISmax[i];
    }
    else
    {
      MD->Ele2IS[i] = MD->ElePrep[i];
      MD->EleIS[i] = MD->EleIS[i] + MD->ElePrep[i]*stepsize; 
    }
    
    MD->EleNetPrep[i] = MD->ElePrep[i] - MD->Ele2IS[i];
  }
}

void calET(realtype t, realtype stepsize, N_Vector VY, void *DS)
{
  int i;
  realtype Delta, Gamma, Es, Ea;
  realtype Rn, G, T, Vel, H, P;
  realtype ET_Value, ET_Remain; 
  realtype *et_Y;
  
  /*
  realtype F_LAI, F_Sita;
  realtype LAI, Sita;
  */
  
  Model_Data MD;
  
  MD = (Model_Data)DS;
  et_Y = NV_DATA_S(VY);
  
  for(i=0; i<MD->NumEle; i++)
  {
    Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn-1], t);
    G = Interpolation(&MD->TSD_G[MD->Ele[i].G-1], t);
    T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp-1], t);
    Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel-1], t);
    H = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity-1], t);
    P = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure-1], t);
    
    Es = 2.53e8*exp(-5.42e3/(T+273));
    Ea = Es*H;
    Gamma = 1e-3*1.013*P/(0.622*2.50036);
    Delta = 4098*Es/(pow(237.3 + T, 2));
    
    /* need to test this */
    MD->EleETP[i] = (1e-3/1440)*(0.408*0.0864*Delta*(Rn - G) + Gamma*900*Vel*(Es - Ea)/(T + 273))/(Delta + Gamma*(1 + 0.34*Vel));
    
    ET_Value = MD->EleETP[i]*stepsize;
    
    // for test cases only;
    /*ET_Value = 0;*/
    
    if (MD->EleIS[i] >= ET_Value)
    {
      MD->EleIS[i] = MD->EleIS[i] - ET_Value;
      MD->EleET[i][0] = ET_Value;
      MD->EleET[i][1] = 0;
      MD->EleET[i][2] = 0;
      MD->EleET[i][3] = 0;
      ET_Remain = 0;
    }
    else
    {
      MD->EleET[i][0] = MD->EleIS[i];
      ET_Remain = ET_Value - MD->EleIS[i];
      MD->EleIS[i] = 0;
      
      if (et_Y[i] >= ET_Remain)
      {
        MD->EleET[i][1] = ET_Remain;
        et_Y[i] = et_Y[i] - ET_Remain;
        MD->EleET[i][2] = 0;
        MD->EleET[i][3] = 0;
        ET_Remain = 0;
      }
      else
      {
        if (et_Y[i] > 0)
        {
          MD->EleET[i][1] = et_Y[i];
          ET_Remain = ET_Remain - et_Y[i];
          et_Y[i] = 0;
        }
        else
        {
          MD->EleET[i][1] = 0;
        }  
        
        /* for shallow groundwater mode */
        if(et_Y[i + 2*MD->NumEle] >= ET_Remain)
        {
          MD->EleET[i][2] = ET_Remain;
          et_Y[i + 2*MD->NumEle] = et_Y[i + 2*MD->NumEle] - ET_Remain/MD->Ele[i].Porosity;
          MD->EleET[i][3] = 0;
        }
        else
        {
          if(et_Y[i + 2*MD->NumEle] >= 0)
          {
            MD->EleET[i][2] = et_Y[i + 2*MD->NumEle];
            et_Y[i + 2*MD->NumEle] = 0;
            MD->EleET[i][3] = ET_Remain - MD->EleET[i][2];
          }
          else
          {
            MD->EleET[i][2] = 0;
            MD->EleET[i][3] = ET_Remain;
          }  
        }  
              
        /*
        LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LAI-1], t);
        Sita = NV_Ith_S(VY, i + MD->NumEle)/(MD->Ele[i].zmax - MD->Ele[i].zmin - NV_Ith_S(VY, i + 2*MD->NumEle));
        
        F_LAI = LAI/8.0;
        
        if (F_LAI > 1.0)
        {
          F_LAI = 1.0;
        }
        else if (F_LAI < 0.0)
        {
          F_LAI = 0.0;
        }
        
        F_Sita =  (Sita - MD->Soil[MD->Ele[i].soil-1].SitaR)/(MD->Soil[MD->Ele[i].soil-1].SitaS - MD->Soil[MD->Ele[i].soil-1].SitaR);
        
        if (F_Sita > 1.0)
        {
          F_Sita = 1.0;
        }
        else if (F_Sita < 0.0)
        {
          F_Sita = 0.0;
        }
        
        if (Sita > MD->Soil[MD->Ele[i].soil-1].SitaR)
        {
          MD->EleET[i][2] = F_LAI*F_Sita*ET_Remain;
         
          NV_Ith_S(VY, i + MD->NumEle) = NV_Ith_S(VY, i + MD->NumEle) - MD->EleET[i][2];

          ET_Remain = ET_Remain - MD->EleET[i][2];
          MD->EleET[i][3] = ET_Remain;
        } 
        */
      } 
    } 
  } 
}
