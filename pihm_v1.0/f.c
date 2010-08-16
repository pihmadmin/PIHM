/*******************************************************************************
 * File        : f.c                                                           *
 * Function    : calculate fluxes, right hand side and construct ODE system    *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               * 
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * This is the model kernel. Idealy, Users can modify this kernel based on     *
 * hydrological conceptual model. Also, flexible constutive relationships can  *
 * apply. This is a multiple scale model in that                               *
 *   1. All equations solved are in integrated form, which make possible to    *
 *      use larger elements;                                                   *
 *   2. Flexible constitutive relationships can be applied;                    *
 *   3. Users can modify this kernel based on appropriate conceptual model.    *
 * We have prove that in finer limit, this model works as good as any other    *
 * tranditional finite element/difference model. However, for larger scale     *
 * application, this model is possible to be applied due to above features.    * 
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

void f(integertype N, realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
{
  int i, j;
  
  realtype Avg_Y_Surf, Dif_Y_Surf;
  realtype Distance;
  realtype Grad_Y_Surf, Avg_Sf;
  
  realtype TotalY_Riv, TotalY_Riv_down;
  realtype Wid, Wid_down;
  realtype Avg_Wid;
  realtype Avg_Y_Riv, Avg_Rough;
  realtype Left_Ele_Y, Left_Ele_YH;
  realtype Right_Ele_Y, Right_Ele_YH;
  realtype Dif_Y_Riv;
  realtype Avg_Y_Sub, Dif_Y_Sub;
  realtype Avg_Ksat, Grad_Y_Sub;
  realtype mp_factor;
  realtype G, GI;
  realtype Cwr, RivPrep;
  realtype temp1, temp2;
  
  realtype Alfa, Beta, CrossA;
  realtype bank_ele;
  realtype AquiferDepth, Deficit, PH;
  
  realtype *Y, *DY;
  Model_Data MD;
  
  Y = NV_DATA_S(CV_Y);
  DY = NV_DATA_S(CV_Ydot);
  MD = (Model_Data) DS;
  
  for(i=0; i<MD->NumEle; i++)
  {
    /*Unconfine Condition */

    if((Y[i + 2*MD->NumEle] >= ((MD->Ele[i].zmax - MD->Ele[i].zmin) + 0.1)))
    {
         Y[i] = Y[i] + MD->Ele[i].Porosity * (Y[i + 2*MD->NumEle] - (MD->Ele[i].zmax - MD->Ele[i].zmin + 0.1));
         Y[i + MD->NumEle] = 0.0;
         Y[i + 2*MD->NumEle] = MD->Ele[i].zmax - MD->Ele[i].zmin + 0.1;   
    }
    
    for(j=0; j<3; j++)
    {
      if(MD->Ele[i].nabr[j] > 0)
      {
        /* groundwater interaction */
        Avg_Y_Sub = (Y[i+2*MD->NumEle] + Y[MD->Ele[i].nabr[j]-1+2*MD->NumEle])/2.0;
        Dif_Y_Sub = (Y[i+2*MD->NumEle] + MD->Ele[i].zmin) - (Y[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle] + MD->Ele[MD->Ele[i].nabr[j]-1].zmin);
        Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x), 2) + pow((MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y), 2));
        Avg_Ksat = (MD->Ele[i].Ksat + MD->Ele[MD->Ele[i].nabr[j] - 1].Ksat)/2.0;
        Grad_Y_Sub = Dif_Y_Sub/Distance;
        
        /* take care of macropore effect */
        if (MD->Soil[(MD->Ele[i].soil-1)].Macropore == 0)
        {
          mp_factor = 1;
        }
        else if (MD->Soil[(MD->Ele[i].soil-1)].Macropore == 1)
        {
          if (Y[i+2*MD->NumEle] > MD->Soil[(MD->Ele[i].soil-1)].base)
          {
            temp1 = pow(10, MD->Soil[(MD->Ele[i].soil-1)].gama*(Y[i+2*MD->NumEle]/MD->Soil[(MD->Ele[i].soil-1)].base - 1));
          }
          else
          {
            temp1 = 1;
          }
          
          if (Y[MD->Ele[i].nabr[j] - 1 + 2*MD->NumEle] > MD->Soil[(MD->Ele[MD->Ele[i].nabr[j]-1].soil-1)].base)
          {
            temp2 = pow(10, MD->Soil[(MD->Ele[MD->Ele[i].nabr[j]-1].soil-1)].gama*(Y[MD->Ele[i].nabr[j]-1+2*MD->NumEle]/MD->Soil[(MD->Ele[MD->Ele[i].nabr[j]-1].soil-1)].base - 1));
          }
          else
          {
            temp2 = 1;
          }  
          
          mp_factor = (temp1 + temp2)/2.0; 
        }
        
        /* groundwater flow modeled by Darcy's law */
        MD->FluxSub[i][j] = mp_factor*Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*MD->Ele[i].edge[j];
        
        /* No Source Check */
        if(Y[i + 2*MD->NumEle] <= 0 && MD->FluxSub[i][j] > 0)
        {
          MD->FluxSub[i][j] = 0;
        }
        if(Y[MD->Ele[i].nabr[j] - 1 + 2*MD->NumEle] <= 0 && MD->FluxSub[i][j] < 0)
        {
          MD->FluxSub[i][j] = 0;
        }
        
        /* Saturation check */
        /*
        if((Y[i + 2*MD->NumEle] >= (MD->Ele[i].zmax - MD->Ele[i].zmin)) && MD->FluxSub[i][j] < 0)
        {
          MD->FluxSub[i][j] = 0;
        }
        if((Y[MD->Ele[i].nabr[j] - 1 + 2*MD->NumEle] >= (MD->Ele[MD->Ele[i].nabr[j] - 1].zmax - MD->Ele[MD->Ele[i].nabr[j] - 1].zmin)) && MD->FluxSub[i][j] <0)
        {
          MD->FluxSub[i][j] = 0;
        }*/
        
        /* Surface Interaction */
        Avg_Y_Surf = (Y[i] + Y[MD->Ele[i].nabr[j]-1])/2.0;
        Dif_Y_Surf = (Y[i] + MD->Ele[i].zmax) - (Y[MD->Ele[i].nabr[j] - 1] + MD->Ele[MD->Ele[i].nabr[j] - 1].zmax);
        Grad_Y_Surf = Dif_Y_Surf/Distance;
        Avg_Sf = (MD->Ele[i].Sf + MD->Ele[MD->Ele[i].nabr[j] - 1].Sf)/2.0;
        Avg_Rough = 0.5*(MD->Ele[i].Rough + MD->Ele[MD->Ele[i].nabr[j] - 1].Rough);
        CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];
        
        /* if surface gradient is not enough to overcome the friction */
        if(fabs(Grad_Y_Surf) <= Avg_Sf)
        {
          MD->FluxSurf[i][j] = 0;
        } 
        else if((Grad_Y_Surf > 0) && (Grad_Y_Surf > Avg_Sf))
        {
          switch(MD->SurfMode) 
          {
            case 1:      
              
              /* Kinematic Wave Approximation constitutive relationship: Manning Equation */
              Alfa = sqrt(Grad_Y_Surf - Avg_Sf)/Avg_Rough;
              Beta = pow(Avg_Y_Surf, 2.0/3.0);
              MD->FluxSurf[i][j] = 60*Alfa*Beta*CrossA;
              break;                   
              
            case 2:
            
              /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */
              Alfa = pow(Avg_Y_Surf, 2.0/3.0)/Avg_Rough;
              Beta = Alfa/sqrt(Grad_Y_Surf - Avg_Sf);
              MD->FluxSurf[i][j] = 60*CrossA*Beta*Grad_Y_Surf;
              break;
              
            default:
              printf("Fatal Error: Surface Overland Mode Type Is Wrong!");
              exit(1);
          }                        
        }
        else if((Grad_Y_Surf < 0) && ((-Grad_Y_Surf) > Avg_Sf))
        {
          switch(MD->SurfMode) 
          {
            case 1:
            
              /* Kinematic Wave Approximation constitutive relationship: Manning Equation */
              Alfa = sqrt(-Grad_Y_Surf - Avg_Sf)/Avg_Rough;
              Beta = pow(Avg_Y_Surf, 2.0/3.0);
              MD->FluxSurf[i][j] = -60*Alfa*Beta*CrossA;
              break;                   
              
            case 2:
              
              /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */
              Alfa = pow(Avg_Y_Surf, 2.0/3.0)/Avg_Rough;
              Beta = Alfa/sqrt(-Grad_Y_Surf - Avg_Sf);
              MD->FluxSurf[i][j] = 60*CrossA*Beta*Grad_Y_Surf;
              break;
              
            default:
              printf("Fatal Error: Surface Overland Mode Type Is Wrong!");
              exit(1);
          }   
        }
        
        /* Source Check */
        if(Y[i] <= 0 && MD->FluxSurf[i][j] > 0)
        {
          MD->FluxSurf[i][j] = 0;
        }
        if(Y[MD->Ele[i].nabr[j] - 1] <= 0 && MD->FluxSurf[i][j] < 0)
        {
          MD->FluxSurf[i][j] = 0;
        }
      }
      else
      {
        /* Handle boundary conditions for elements. No flow (natural) boundary condition is default */
        if(MD->Ele[i].BC == 0)
        {
          MD->FluxSurf[i][j] = 0;
          MD->FluxSub[i][j] = 0;
        }
        else
        {
          MD->FluxSurf[i][j] = 0;
          
          /* this part of code has not been tested! */
          if(MD->Ele[i].BC > 0)         /* Dirichlet BC */
          {
            Avg_Y_Sub = (Y[i+2*MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC)-1], t) - MD->Ele[i].zmin))/2;
            Dif_Y_Sub = (Y[i+2*MD->NumEle] + MD->Ele[i].zmin) - 
                         Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC)-1], t);            
            Distance = sqrt(pow(MD->Ele[i].edge[0]*MD->Ele[i].edge[1]*MD->Ele[i].edge[2]/(4*MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j]/2, 2));
            Avg_Ksat = MD->Ele[i].Ksat;
            Grad_Y_Sub = Dif_Y_Sub/Distance;
        
            MD->FluxSub[i][j] = Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*MD->Ele[i].edge[j];
          }
          else                          /* Nuemann BC */
          {
            MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC)-1+MD->Num1BC], t);
          }
        }
        
        /* Source check */
        if(Y[i + 2*MD->NumEle] <= 0 && MD->FluxSub[i][j] > 0)
        {
          MD->FluxSub[i][j] = 0;
        }  
      }
    }
  }
  
  /* initialize river flux */
  for(i=0; i<MD->NumRiv; i++)
  {
    for(j=0; j<6; j++)
    {
      MD->FluxRiv[i][j] = 0;
    }  
  } 
   
  for(i=0; i<MD->NumRiv; i++)
  {
    /* Note: the ordering of river segment in input file has to be: from UP to DOWN */   
    TotalY_Riv = Y[i + 3*MD->NumEle] + MD->Riv[i].zmin;
    Wid = MD->Riv_Shape[MD->Riv[i].shape - 1].width;
    
    if(MD->Riv[i].down > 0)
    {
      TotalY_Riv_down = Y[MD->Riv[i].down - 1 + 3*MD->NumEle] + MD->Riv[MD->Riv[i].down - 1].zmin;
      Wid_down = MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].width;
      Avg_Wid = (Wid + Wid_down)/2.0;
      Avg_Y_Riv = (Y[i+3*MD->NumEle] + Y[MD->Riv[i].down - 1 + 3*MD->NumEle])/2.0;
      Avg_Rough = (MD->Riv_Mat[MD->Riv[i].material - 1].Rough + MD->Riv_Mat[MD->Riv[MD->Riv[i].down - 1].material-1].Rough)/2.0;
      Distance = sqrt(pow(MD->Riv[i].x - MD->Riv[MD->Riv[i].down - 1].x, 2) + pow(MD->Riv[i].y - MD->Riv[MD->Riv[i].down - 1].y, 2));
      Dif_Y_Riv = (TotalY_Riv - TotalY_Riv_down)/Distance; 
      Avg_Sf = (MD->Riv_Mat[MD->Riv[i].material - 1].Sf + MD->Riv_Mat[MD->Riv[MD->Riv[i].down - 1].material-1].Sf)/2.0;   
      CrossA = Avg_Y_Riv*Avg_Wid;
      
      if(fabs(Dif_Y_Riv) <= Avg_Sf)
      {
        MD->FluxRiv[i][1] = 0;
      } 
      else if((Dif_Y_Riv > 0) && (Dif_Y_Riv > Avg_Sf))
      {
        switch(MD->RivMode) 
        {
          case 1:  
          
          /* Kinematic Wave Approximation constitutive relationship: Manning Equation */    
          Alfa = sqrt(Dif_Y_Riv - Avg_Sf)/(Avg_Rough*pow((Avg_Wid + 2*Avg_Y_Riv), 2.0/3.0));
          Beta = 5.0/3.0;
          MD->FluxRiv[i][1] = 60*Alfa*pow(CrossA, Beta);
          break;                   
              
          case 2:
          
          /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */    
          Alfa = pow(Avg_Y_Riv, 2.0/3.0)/Avg_Rough;
          Beta = Alfa/sqrt(Dif_Y_Riv - Avg_Sf);
          MD->FluxRiv[i][1] = 60*CrossA*Beta*Dif_Y_Riv;
          break;
              
          default:
            printf("Fatal Error: River Routing Mode Type Is Wrong!");
            exit(1);
        }
      }
      else if((Dif_Y_Riv < 0) && (-Dif_Y_Riv > Avg_Sf))
      {
        switch(MD->RivMode) 
        {
          case 1:
              
          /* Kinematic Wave Approximation constitutive relationship: Manning Equation */    
          Alfa = sqrt(-Dif_Y_Riv - Avg_Sf)/(Avg_Rough*pow((Avg_Wid + 2*Avg_Y_Riv), 2.0/3.0));
          Beta = 5.0/3.0;
          MD->FluxRiv[i][1] = -60*Alfa*pow(CrossA, Beta);
          break;                     
              
          case 2:
          
          /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */     
          Alfa = pow(Avg_Y_Riv, 2.0/3.0)/Avg_Rough;
          Beta = Alfa/sqrt(-Dif_Y_Riv - Avg_Sf);
          MD->FluxRiv[i][1] = 60*CrossA*Beta*Dif_Y_Riv;
          break;
              
          default:
            printf("Fatal Error: River Routing Mode Type Is Wrong!");
            exit(1);
        }
      }  

      /* Source check */
      if(Y[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] > 0)
      {
        MD->FluxRiv[i][1] = 0;
      }
      else if(Y[MD->Riv[i].down - 1 + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] < 0)
      {
        MD->FluxRiv[i][1] = 0;
      }
      
      /* accumulate to get in-flow for down segments */
      MD->FluxRiv[MD->Riv[i].down - 1][0] = MD->FluxRiv[MD->Riv[i].down - 1][0] + MD->FluxRiv[i][1];
    }
    else
    {
      switch(MD->Riv[i].down)
      {
        case -1:  
         
          /* Dirichlet boundary condition */        
          TotalY_Riv_down = Interpolation(&MD->TSD_Riv[(MD->Riv[i].BC)-1], t) + MD->Node[MD->Riv[i].ToNode-1].zmin + MD->Riv_Shape[MD->Riv[i].shape-1].bed;
          Distance = sqrt(pow(MD->Riv[i].x - MD->Node[MD->Riv[i].ToNode-1].x, 2) + pow(MD->Riv[i].y - MD->Node[MD->Riv[i].ToNode-1].y, 2));
          Dif_Y_Riv = (TotalY_Riv - TotalY_Riv_down)/Distance; 
          Avg_Sf = MD->Riv_Mat[MD->Riv[i].material - 1].Sf;
          Avg_Rough = MD->Riv_Mat[MD->Riv[i].material-1].Rough;
          Avg_Y_Riv = Y[i + 3*MD->NumEle];
          Avg_Wid = Wid;
          CrossA = Wid*Avg_Y_Riv;
          
          if(fabs(Dif_Y_Riv) <= Avg_Sf)
          {
            MD->FluxRiv[i][1] = 0;
          } 
          else if((Dif_Y_Riv > 0) && (Dif_Y_Riv > Avg_Sf))
          {
            switch(MD->RivMode) 
            {
              case 1:  
                
                /* Kinematic Wave Approximation constitutive relationship: Manning Equation */    
                Alfa = sqrt(Dif_Y_Riv - Avg_Sf)/(Avg_Rough*pow((Avg_Wid + 2*Avg_Y_Riv), 2.0/3.0));
                Beta = 5.0/3.0;
                MD->FluxRiv[i][1] = 60*Alfa*pow(CrossA, Beta);
                break;                
              
              case 2:     
              
                /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */
                Alfa = pow(Avg_Y_Riv, 2.0/3.0)/Avg_Rough;
                Beta = Alfa/sqrt(Dif_Y_Riv - Avg_Sf);
                MD->FluxRiv[i][1] = 60*CrossA*Beta*Dif_Y_Riv;
                break;
              
              default:
                printf("Fatal Error: River Routing Mode Type Is Wrong!");
                exit(1);
            }
          }
          else if((Dif_Y_Riv < 0) && (-Dif_Y_Riv > Avg_Sf))
          {
            switch(MD->RivMode) 
            {
              case 1: 
              
                /* Kinematic Wave Approximation constitutive relationship: Manning Equation */    
                Alfa = sqrt(-Dif_Y_Riv - Avg_Sf)/(Avg_Rough*pow((Avg_Wid + 2*Avg_Y_Riv), 2.0/3.0));
                Beta = 5.0/3.0;
                MD->FluxRiv[i][1] = -60*Alfa*pow(CrossA, Beta);
                break;                      
              
              case 2:
              
                /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */
                Alfa = pow(Avg_Y_Riv, 2.0/3.0)/Avg_Rough;
                Beta = Alfa/sqrt(-Dif_Y_Riv - Avg_Sf);
                MD->FluxRiv[i][1] = 60*CrossA*Beta*Dif_Y_Riv;
                break;
              
              default:
                printf("Fatal Error: River Routing Mode Type Is Wrong!");
                exit(1);
            }
          }
          break;
          
        case -2:
          
          /* Neumann boundary condition */
          MD->FluxRiv[i][1] = Interpolation(&MD->TSD_Riv[MD->Riv[i].BC-1], t);
          break;
          
        case -3:    
          
          /* zero-depth-gradient boundary conditions */
          Distance = sqrt(pow(MD->Riv[i].x - MD->Node[MD->Riv[i].ToNode-1].x, 2) + pow(MD->Riv[i].y - MD->Node[MD->Riv[i].ToNode-1].y, 2));
          Dif_Y_Riv = (MD->Riv[i].zmin - (MD->Node[MD->Riv[i].ToNode-1].zmin + MD->Riv_Shape[MD->Riv[i].shape-1].bed))/Distance;
          Avg_Rough = MD->Riv_Mat[MD->Riv[i].material-1].Rough;
          Avg_Y_Riv = Y[i + 3*MD->NumEle];
          Avg_Wid = Wid;
          CrossA = Wid*Avg_Y_Riv;
          MD->FluxRiv[i][1] = 60*Avg_Wid*pow(Avg_Y_Riv, 5.0/3.0)*sqrt(Dif_Y_Riv)/Avg_Rough;
          break;
          
        case -4:
          
          /* Critical Depth boundary conditions */
          CrossA = Wid*Y[i + 3*MD->NumEle];
          MD->FluxRiv[i][1] = 60*CrossA*sqrt(9.81*Y[i + 3*MD->NumEle]);
          break;
          
        default:
          printf("Fatal Error: River Routing Boundary Condition Type Is Wrong!");
          exit(1); 
      }     
      
      /* Source check */
      if(Y[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] > 0)
      {
        MD->FluxRiv[i][1] = 0;
      }
      else if(Y[MD->Riv[i].down - 1 + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] < 0)
      {
        MD->FluxRiv[i][1] = 0;
      }  
      
      /* need some work to take care of multiple output Q */
      MD->Q = MD->FluxRiv[i][1];
          
    }
    
    /* Interaction between surface flow and channel */
    if (MD->Riv[i].LeftEle > 0)
    {
      Left_Ele_Y = Y[MD->Riv[i].LeftEle - 1];
      Left_Ele_YH = Y[MD->Riv[i].LeftEle - 1] + MD->Ele[MD->Riv[i].LeftEle - 1].zmax;
      Distance = sqrt(pow(MD->Riv[i].x - MD->Ele[MD->Riv[i].LeftEle - 1].x, 2) + pow(MD->Riv[i].y - MD->Ele[MD->Riv[i].LeftEle - 1].y, 2));
      Cwr = MD->Riv_Mat[MD->Riv[i].material-1].Cwr;
      
      if (MD->Riv[i].zmax < MD->Ele[MD->Riv[i].LeftEle - 1].zmax) { bank_ele = MD->Ele[MD->Riv[i].LeftEle - 1].zmax;} 
      else { bank_ele = MD->Riv[i].zmax;}
      
      if (TotalY_Riv > Left_Ele_YH)
      {
        if (Left_Ele_YH > bank_ele)
        {
          MD->FluxRiv[i][2] = Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(TotalY_Riv - Left_Ele_YH)*(TotalY_Riv - bank_ele)/3.0;
        }
        else
        {
          MD->FluxRiv[i][2] = Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(TotalY_Riv - bank_ele)*(TotalY_Riv - bank_ele)/3.0;
        }
      }
      else
      {
        if (TotalY_Riv > bank_ele)
        {
          MD->FluxRiv[i][2] = -Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(Left_Ele_YH - TotalY_Riv)*(Left_Ele_YH - bank_ele)/3.0;
        }
        else
        {
          /* This is basicaly a lump representation */
          MD->FluxRiv[i][2] = -Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(Left_Ele_YH - bank_ele)*(Left_Ele_YH - bank_ele)/3.0;
        }
      }
      
      /* source check */
      if(Y[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][2] > 0)
      {
        MD->FluxRiv[i][2] = 0;
      }
    
      if(Y[MD->Riv[i].LeftEle - 1] <= 0 && MD->FluxRiv[i][2] < 0)
      {
        MD->FluxRiv[i][2] = 0;
      }
      
      /* replace overland flux item */
      for(j=0; j < 3; j++)
      {
        /* this may cause trouble when there are 2 boundary side and one is not neutral.*/
        if(MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j] == MD->Riv[i].RightEle)
        {
          MD->FluxSurf[MD->Riv[i].LeftEle - 1][j] = -MD->FluxRiv[i][2];
        }
      }      
    }
    
    if (MD->Riv[i].RightEle > 0)
    {
      Right_Ele_Y = Y[MD->Riv[i].RightEle - 1];
      Right_Ele_YH = Y[MD->Riv[i].RightEle - 1] + MD->Ele[MD->Riv[i].RightEle - 1].zmax;
      Distance = sqrt(pow(MD->Riv[i].x - MD->Ele[MD->Riv[i].RightEle - 1].x, 2) + pow(MD->Riv[i].y - MD->Ele[MD->Riv[i].RightEle - 1].y, 2));
      Cwr = MD->Riv_Mat[MD->Riv[i].material-1].Cwr;
      
      if (MD->Riv[i].zmax < MD->Ele[MD->Riv[i].RightEle - 1].zmax) { bank_ele = MD->Ele[MD->Riv[i].RightEle - 1].zmax;} 
      else { bank_ele = MD->Riv[i].zmax;}
      
      if (TotalY_Riv > Right_Ele_YH)
      {
        if (Right_Ele_YH > bank_ele)
        {
          MD->FluxRiv[i][3] = Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(TotalY_Riv - Right_Ele_YH)*(TotalY_Riv - bank_ele)/3.0;
        }
        else
        {
          MD->FluxRiv[i][3] = Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(TotalY_Riv - bank_ele)*(TotalY_Riv - bank_ele)/3.0;
        }
      }
      else
      {
        if (TotalY_Riv > bank_ele)
        {
          MD->FluxRiv[i][3] = -Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(Right_Ele_YH - TotalY_Riv)*(Right_Ele_YH - bank_ele)/3.0;
        }
        else
        {
          MD->FluxRiv[i][3] = -Cwr*60*2.0*sqrt(2*9.81)*MD->Riv[i].Length*sqrt(Right_Ele_YH - bank_ele)*(Right_Ele_YH - bank_ele)/3.0;
        }
      }
      
      /* source check */
      if(Y[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][3] > 0)
      {
        MD->FluxRiv[i][3] = 0;
      }
    
      if(Y[MD->Riv[i].RightEle - 1] <= 0 && MD->FluxRiv[i][3] < 0)
      {
        MD->FluxRiv[i][3] = 0;
      }
      
      /* replace overland flux item */
      for(j=0; j < 3; j++)
      {
        if(MD->Ele[MD->Riv[i].RightEle - 1].nabr[j] == MD->Riv[i].LeftEle)
        {
          MD->FluxSurf[MD->Riv[i].RightEle - 1][j] = -MD->FluxRiv[i][3];
          break;
        }
      }
    }
           
    /* groundwater interaction */
    if (MD->Riv[i].LeftEle > 0)
    {
      /* Left Neighbor Groundwater Head */
      Left_Ele_Y = Y[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle];
      Left_Ele_YH = Y[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle - 1].zmin;
      Distance = sqrt(pow(MD->Riv[i].x - MD->Ele[MD->Riv[i].LeftEle - 1].x, 2) + pow(MD->Riv[i].y - MD->Ele[MD->Riv[i].LeftEle - 1].y, 2));
    
      if (MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].Macropore == 0)
      {
        mp_factor = 1;
      }
      else if (MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].Macropore == 1)
      {
        if (Y[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] > MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].base)
        {
          mp_factor = pow(10, MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].gama*(Y[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle]/MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].base - 1));
        }
        else
        {
          mp_factor = 1;
        }
      }
    
      MD->FluxRiv[i][4] = mp_factor*MD->Riv[i].Length*(0.5*Wid+Y[i+3*MD->NumEle]) * MD->Ele[MD->Riv[i].LeftEle - 1].Ksat * (TotalY_Riv - Left_Ele_YH)/Distance;
    
      /* source check */
      if(Y[i+3*MD->NumEle] <= 0 && MD->FluxRiv[i][4] > 0)
      {
        MD->FluxRiv[i][4] = 0;
      }
      if(Y[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] <= 0 && MD->FluxRiv[i][4] < 0)
      {
        MD->FluxRiv[i][4] = 0;
      }
    }
    
    if (MD->Riv[i].RightEle > 0)  
    {
      Right_Ele_Y = Y[MD->Riv[i].RightEle - 1 + 2*MD->NumEle];
      Right_Ele_YH = Y[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].RightEle - 1].zmin;
      Distance = sqrt(pow(MD->Riv[i].x - MD->Ele[MD->Riv[i].RightEle - 1].x, 2) + pow(MD->Riv[i].y - MD->Ele[MD->Riv[i].RightEle - 1].y, 2));
    
      if (MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].Macropore == 0)
      {
        mp_factor = 1;
      }
      else if (MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].Macropore == 1)
      {
        if (Y[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] > MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].base)
        {
          mp_factor = pow(10, MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].gama*(Y[MD->Riv[i].RightEle - 1 + 2*MD->NumEle]/MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].base - 1));
        }
        else
        {
          mp_factor = 1;
        }
      }
      
      MD->FluxRiv[i][5] = mp_factor*MD->Riv[i].Length*(0.5*Wid+Y[i+3*MD->NumEle]) * MD->Ele[MD->Riv[i].RightEle - 1].Ksat * (TotalY_Riv - Right_Ele_YH)/Distance;
    
      /* source check */
      if(Y[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][5] > 0)
      {
        MD->FluxRiv[i][5] = 0;
      }
      if(Y[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] <= 0 && MD->FluxRiv[i][5] < 0)
      {
        MD->FluxRiv[i][5] = 0;
      }
    }  
  }
  
  /* Calculate DY */ 
  switch(MD->UnsatMode) {
  
  case 1: 
   
    /* Shallow Groundwater assumption applied here*/
    for(i=0; i<MD->NumEle; i++)
    {
      AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
      
      if(Y[i+2*MD->NumEle] >= AquiferDepth)
      {
        Deficit = 0;
        MD->EleVic[i] = 0;
      }
      else
      {
        Deficit = AquiferDepth - Y[i+2*MD->NumEle];
        MD->EleVic[i] = Interpolation(&MD->TSD_Inc[MD->Soil[(MD->Ele[i].soil-1)].Inf-1], t);
      } 

      /* handle Prep and Infiltration first */
      if(Y[i+MD->NumEle] <= Deficit)
      {
        if(Y[i] > 0)
        {
          DY[i] = MD->EleNetPrep[i] - MD->EleVic[i];
          DY[i+2*MD->NumEle] = MD->EleVic[i];
        } 
        else if(MD->EleNetPrep[i] > MD->EleVic[i] && Y[i] <= 0)
        {
          DY[i] = MD->EleNetPrep[i] - MD->EleVic[i];
          DY[i+2*MD->NumEle] = MD->EleVic[i];
        }
        else if(MD->EleNetPrep[i] < MD->EleVic[i] && Y[i] <= 0 && MD->EleNetPrep[i] > 0)
        {
          DY[i] = 0;
          DY[i+2*MD->NumEle] = MD->EleNetPrep[i];
        }
        else     
        {
          DY[i] = 0;
          DY[i+2*MD->NumEle] = MD->EleNetPrep[i];
        }
      }
      else
      {
        /* The reason of this happening is not clear, therefore the treatment is not secure
           Fortranately, this will not happen in most cases. */
        DY[i] = MD->EleNetPrep[i];
        DY[i+2*MD->NumEle] = 0;
 
        if (Deficit > 0)
        {
          Y[i+MD->NumEle] = AquiferDepth - Y[i+2*MD->NumEle];
        }
        else
        {
          Y[i+2*MD->NumEle] = AquiferDepth;
          Y[i+MD->NumEle] = 0;
        }  
      }
      
      /* handle surface flux then */
      for(j=0; j<3; j++)
      {
        DY[i] =  DY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
      }  
      
      /*
      if(Y[i] <= 0 && DY[i] < 0)
      {
        DY[i] = 0;
        DY[i+2*MD->NumEle] = MD->EleNetPrep[i];
      }*/
      
      for(j=0; j<3; j++)
      {
        DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle] - MD->FluxSub[i][j]/MD->Ele[i].area;
      }
    }
  
    for(i=0; i<MD->NumRiv; i++)
    {
      if (MD->Riv[i].LeftEle > 0 && MD->Riv[i].RightEle > 0)
      {
        RivPrep = Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].LeftEle-1].prep-1], t) + Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].RightEle-1].prep-1], t);
        RivPrep = RivPrep/2;
      }  
      else if (MD->Riv[i].LeftEle > 0 && MD->Riv[i].RightEle == 0)
      {
        RivPrep = Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].LeftEle-1].prep-1], t);
      }
      else if (MD->Riv[i].LeftEle == 0 && MD->Riv[i].RightEle > 0)  
      {
        RivPrep = Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].RightEle-1].prep-1], t);
      }
      else
      {  
        /* for test case 3 only */
        /*RivPrep = Interpolation(&MD->TSD_Prep[0], t);*/
        RivPrep = 0;
      }
      
      DY[i+3*MD->NumEle] = RivPrep + MD->FluxRiv[i][0] - MD->FluxRiv[i][1];
      DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle] - MD->FluxRiv[i][2] - MD->FluxRiv[i][3];
      DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle] - MD->FluxRiv[i][4] - MD->FluxRiv[i][5];
      DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle]/(MD->Riv[i].Length*Wid);
    
      if (MD->Riv[i].LeftEle > 0)
      {
        DY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] = DY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]
                                + MD->FluxRiv[i][4]/MD->Ele[MD->Riv[i].LeftEle-1].area;
      }
      if (MD->Riv[i].RightEle > 0)
      {                          
        DY[MD->Riv[i].RightEle-1 + 2*MD->NumEle] = DY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]
                                 + MD->FluxRiv[i][5]/MD->Ele[MD->Riv[i].RightEle-1].area;
      }
    }
  
    for(i=0; i<MD->NumEle; i++)
    { 
      AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
      Deficit = AquiferDepth - Y[i+2*MD->NumEle];
      G = 1e-4 + MD->Ele[i].Porosity*(1-pow( 1+pow(MD->Ele[i].Alpha*Deficit , MD->Ele[i].Beta),-(MD->Ele[i].Beta+1)/MD->Ele[i].Beta));
      GI = -pow( 1+pow(MD->Ele[i].Alpha*Deficit , MD->Ele[i].Beta), -(MD->Ele[i].Beta+1)/MD->Ele[i].Beta);
    
      /* for test case 1 only */
      /*
      G = MD->Ele[i].Porosity;
      GI = 0; 
      */
      
      DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle]/G;
      DY[i+MD->NumEle] = GI*DY[i+2*MD->NumEle];
      
      // if there is a source (well) in it
      if(MD->Ele[i].source > 0)
      {
        DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle] - Interpolation(&MD->TSD_Source[MD->Ele[i].source - 1], t)/(MD->Ele[i].Porosity*MD->Ele[i].area);
      }
    
      // check if it is out of bound
      if(Y[i+MD->NumEle]>Deficit && DY[i+MD->NumEle]>0)
      {
        DY[i+MD->NumEle] = 0;
      }
      if(Y[i+MD->NumEle]<0 && DY[i+MD->NumEle]<0)
      {
        DY[i+MD->NumEle] = 0;
      }
    
      if(Y[i+2*MD->NumEle]>AquiferDepth && DY[i+2*MD->NumEle]>0)
      {
        DY[i+2*MD->NumEle] = 0;
      }
      if(Y[i+2*MD->NumEle]<0 && DY[i+2*MD->NumEle]<0)
      {
        DY[i+2*MD->NumEle] = 0;
      }
    }  
  
    break;
    
  case 2:
    
    for(i=0; i<MD->NumEle; i++)
    {
      /* Uncomment if there is no interception storage */
      /* MD->EleNetPrep[i] = Interpolation(&MD->TSD_Prep[MD->Ele[i].prep-1], t); */
      MD->EleVic[i] = Interpolation(&MD->TSD_Inc[MD->Soil[(MD->Ele[i].soil-1)].Inf-1], t);
  
      AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
      Deficit = AquiferDepth - Y[i+2*MD->NumEle];

      if(Y[i+MD->NumEle] < Deficit)
      {
        if(Y[i] > 0)
        {
          DY[i] = MD->EleNetPrep[i] - MD->EleVic[i];
          DY[i+MD->NumEle] = MD->EleVic[i];
        }
        else if(MD->EleNetPrep[i] > MD->EleVic[i] && Y[i] <= 0)
        {
          DY[i] = MD->EleNetPrep[i] - MD->EleVic[i];
          DY[i+MD->NumEle] = MD->EleVic[i];
        }
        else if(MD->EleNetPrep[i] < MD->EleVic[i] && Y[i] <= 0 && MD->EleNetPrep[i] > 0)
        {
          DY[i] = 0;
          DY[i+MD->NumEle] = MD->EleNetPrep[i];
        }
        else
        {
          DY[i] = 0;
          DY[i+MD->NumEle] = MD->EleNetPrep[i];
        }  
      }
      else
      {
        DY[i] = MD->EleNetPrep[i];
        DY[i+MD->NumEle] = 0;
      }
      
      for(j=0; j<3; j++)
      {
        DY[i] =  DY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
      }
      
      if(Y[i] <= 0 && DY[i] < 0)
      {
        DY[i] = 0;
        DY[i+MD->NumEle] = MD->EleNetPrep[i];
      }
    
      PH = 1 - exp(-MD->Ele[i].Ksat*Deficit);
     
      MD->Recharge[i] = MD->Ele[i].Ksat*(PH - MD->Ele[i].Alpha*Y[i+MD->NumEle])/(1e-7 + MD->Ele[i].Alpha*Deficit - PH);
    
      if(Y[i+MD->NumEle]<0 && MD->Recharge[i] < 0)
      {
        MD->Recharge[i] = 0;
      }
      if(Y[i+2*MD->NumEle]<0 && MD->Recharge[i] > 0)
      {
        MD->Recharge[i] = 0;
      }
    
      DY[i+MD->NumEle] = DY[i+MD->NumEle] + MD->Recharge[i];
      DY[i+MD->NumEle] = DY[i+MD->NumEle]/MD->Ele[i].Porosity;
      
      /* Source check */
      if(Y[i+MD->NumEle]>Deficit && DY[i+MD->NumEle]>0)
      {
        DY[i+MD->NumEle] = 0;
      }
      if(Y[i+MD->NumEle]<0 && DY[i+MD->NumEle]<0)
      {
        DY[i+MD->NumEle] = 0;
      }
    
      DY[i+2*MD->NumEle] = -MD->Recharge[i];
    
      for(j=0; j<3; j++)
      {
        DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle] - MD->FluxSub[i][j]/MD->Ele[i].area;
      }
    }
  
    for(i=0; i<MD->NumRiv; i++)
    {
      if (MD->Riv[i].LeftEle > 0 && MD->Riv[i].RightEle > 0)
      {
        RivPrep = Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].LeftEle-1].prep-1], t) + Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].RightEle-1].prep-1], t);
        RivPrep = RivPrep/2;
      }  
      else if (MD->Riv[i].LeftEle > 0 && MD->Riv[i].RightEle == 0)
      {
        RivPrep = Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].LeftEle-1].prep-1], t);
      }
      else if (MD->Riv[i].LeftEle == 0 && MD->Riv[i].RightEle > 0)  
      {
        RivPrep = Interpolation(&MD->TSD_Prep[MD->Ele[MD->Riv[i].RightEle-1].prep-1], t);
      }
      else
      {  
        /* for test case 3 only */
        /*RivPrep = Interpolation(&MD->TSD_Prep[0], t);*/
        RivPrep = 0;
      }
      
      DY[i+3*MD->NumEle] = RivPrep + MD->FluxRiv[i][0] - MD->FluxRiv[i][1];
      DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle] - MD->FluxRiv[i][2] - MD->FluxRiv[i][3];
      DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle] - MD->FluxRiv[i][4] - MD->FluxRiv[i][5];
      DY[i+3*MD->NumEle] = DY[i+3*MD->NumEle]/(MD->Riv[i].Length*MD->Riv_Shape[MD->Riv[i].shape-1].width);
    
      if (MD->Riv[i].LeftEle > 0)
      {
        DY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] = DY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] + MD->FluxRiv[i][4]/MD->Ele[MD->Riv[i].LeftEle-1].area;
      }
      
      if (MD->Riv[i].RightEle > 0)
      {                         
        DY[MD->Riv[i].RightEle-1 + 2*MD->NumEle] = DY[MD->Riv[i].RightEle-1 + 2*MD->NumEle] + MD->FluxRiv[i][5]/MD->Ele[MD->Riv[i].RightEle-1].area;
      }                           
    }
  
    for(i=0; i<MD->NumEle; i++)
    {
      
      /* if there is a source (well) in it */
      if(MD->Ele[i].source > 0)
      {
        DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle] - Interpolation(&MD->TSD_Source[MD->Ele[i].source - 1], t)/MD->Ele[i].area;
      }
      
      DY[i+2*MD->NumEle] = DY[i+2*MD->NumEle]/MD->Ele[i].Porosity;
    
      if(Y[i+2*MD->NumEle]>AquiferDepth && DY[i+2*MD->NumEle]>0)
      {
        DY[i+2*MD->NumEle] = 0;
      }
      if(Y[i+2*MD->NumEle]<0 && DY[i+2*MD->NumEle]<0)
      {
        DY[i+2*MD->NumEle] = 0;
      }
    }  
  
    break;
    
  default:
  
    printf("Fatal Error: Unsaturated Layer Mode Type Is Wrong!");
    exit(1);
  
  }
}


realtype Interpolation(TSD *Data, realtype t)
{
  int i, success;
  realtype result;
  
  i=0;
  success = 0;
  
  while(i<Data->length && t>Data->TS[i][0])
  {
   i++;
  }
  
  if(i==0)
  {
    /* t is smaller than the 1st node */
    result = Data->TS[i][1];
  }
  else if(i >= Data->length)
  {
    result = Data->TS[i-1][1];
  }
  else
  {
    result = ((Data->TS[i][0]-t)*Data->TS[i-1][1] + (t-Data->TS[i-1][0])*Data->TS[i][1])/(Data->TS[i][0]-Data->TS[i-1][0]);
    success = 1;
  }
  
  if(success == 0)
  {
    /*
    printf("\nWarning:  Exterpolation is used ...\n");
    */
  }
  
  return result;
}

