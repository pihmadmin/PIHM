/*******************************************************************************
 * File        : f.c                                                           *
 * Function    : calculate fluxes, right hand side and construct ODE system    *
 * Programmers : Yizhong Qu   @ Pennsylvania State Univeristy                  *
 *               Mukesh Kumar @ Pennsylvania State Univeristy                  *
 *               Gopal Bhatt  @ Pennsylvania State Univeristy                  *
 * Version     : 2.0 (July 10, 2007)                                           *
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

//! @file f.c calculate fluxes, right hand side and construct ODE system

/*    C Header Files    */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*    SUNDIALS Header Files    */
#include "nvector_serial.h"
#include "sundials_types.h"

/*    PIHM Header Files    */
#include "pihm.h"
#include "calib.h"

#define EPSILON 0.05
#define THRESH  0.00


#define GRAV  73231257600.0	/**< Gravitational Accelaration in m/day^2 */

#define ABS_TOL    1E-4		/**< Absolute tolerance as defined in control data */

realtype Vic_CALIB;
realtype rivK_CALIB;
realtype Kh_CALIB;
realtype Rec_CALIB;
realtype et2_CALIB;
realtype et1_CALIB;
realtype sat_THRESH;
realtype mp_MULTFH;
realtype mp_MULTFV;
realtype mpArea_CALIB;
realtype ovl_THRESH_H;
realtype ovl_THRESH_V;
realtype rzd_CALIB;


/*    Function Declarations    */
realtype Interpolation(TSD *Data, realtype t);
realtype returnVal(realtype rArea, realtype rPerem, realtype eqWid, realtype ap_Bool);
realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool);
void OverlandFlow(realtype **flux, int loci, int locj, int surfmode, realtype avg_y, realtype grad_y, realtype avg_sf, realtype alfa, realtype beta, realtype crossA, realtype avg_rough, int eletypeBool, realtype avg_perem);
void OLflowFromEleToRiv(realtype sideEle_y,realtype sideEle_zmax,realtype rivX,realtype sideEleX,realtype rivY,realtype sideEleY,realtype cwr,realtype rivZmax,realtype loc_yriver,realtype **fluxriv,int loc_i,int loc_j,realtype length);
void GWflowFromEleToRiv(realtype sideEle_y,realtype sideEle_zmax, realtype sideEle_zmin,realtype rivX,realtype sideEleX,realtype rivY,realtype sideEleY,int loc_McPore,realtype loc_yriver,realtype loc_totyriver,realtype **fluxriv,int loc_i,int loc_j,realtype length, realtype loc_base,realtype loc_gama, realtype loc_perem,realtype loc_ksat,realtype ele_Thresh);







int f(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
//! Function calculates fluxes, right hand side leading to construction of ODE system
/*! \param t is the time of simulation
    \param CV_Y	is state variable vector
    \param CV_Ydot	is vector of rate of change of state variables
    \param DS is pointer to model data structure
*/
{
    int i, j;

    realtype Delta, Gamma;
    realtype Rn, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a,r_s,alpha_r,f_r,eta_s,beta_s,Rmax;


    realtype Avg_Y_Surf, Dif_Y_Surf;
    realtype Distance;
    realtype Grad_Y_Surf, Avg_Sf;

    realtype TotalY_Riv, TotalY_Riv_down;
    realtype Perem, Perem_down, loc_perem;
    realtype Avg_Perem;
    realtype Avg_Y_Riv, Avg_Rough;
    realtype Left_Ele_Y, Left_Ele_YH;
    realtype Right_Ele_Y, Right_Ele_YH;
    realtype Dif_Y_Riv;
    realtype Avg_Y_Sub, Dif_Y_Sub;
    realtype Avg_Ksat, Grad_Y_Sub;
    realtype mp_factor,m_factor;
    realtype G, GI;
    realtype Cwr, RivPrep;
    realtype temp1, temp2;

    realtype Alfa, Beta, CrossA;
    realtype bank_ele;
    realtype AquiferDepth, Deficit, PH,elemSatn,eleSatn;
    realtype loc_bcEle, Avg_BedDepth;

    realtype *Y, *DY,*DummyY,*DummyDY;
    Model_Data MD;

    Vic_CALIB=setVic_CALIB();
    rivK_CALIB=setrivK_CALIB();
    Kh_CALIB=setKh_CALIB();
    Rec_CALIB=setRec_CALIB();
    et2_CALIB=setet2_CALIB();
    et1_CALIB=setet1_CALIB();
    sat_THRESH=setsat_THRESH();
    mp_MULTFH=setmp_MULTFH();
    mp_MULTFV=setmp_MULTFV();
    mpArea_CALIB=setmpArea_CALIB();
    ovl_THRESH_H=setovl_THRESH_H();
    ovl_THRESH_V=setovl_THRESH_V();
    rzd_CALIB=setrzd_CALIB();


    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
    MD = (Model_Data) DS;

    DummyY=(realtype *)malloc((3*MD->NumEle+MD->NumRiv)*sizeof(realtype));
    DummyDY=(realtype *)malloc((3*MD->NumEle+MD->NumRiv)*sizeof(realtype));


    /* Lateral Flux Calculation Follows */
    for(i=0; i<3*MD->NumEle+MD->NumRiv; i++)
    {

        //printf("%.35lf %d %lf\n",t, i, Y[i]);

         DummyDY[i]=0.0;
           if(Y[i]<=ABS_TOL)
           {
              DummyY[i]=0;
         }
           else
           {
                DummyY[i]=Y[i];
         }
    }
    for(i=MD->NumEle; i<2*MD->NumEle; i++)
    {
          //      if((i>=MD->NumEle)&&(i<2*MD->NumEle)){
        if(Y[i]>(MD->Ele[i-MD->NumEle].zmax-MD->Ele[i-MD->NumEle].zmin))
        {
            DummyY[i]=MD->Ele[i-MD->NumEle].zmax - MD->Ele[i-MD->NumEle].zmin;
        }
         if(Y[i+MD->NumEle]>(MD->Ele[i-MD->NumEle].zmax-MD->Ele[i-MD->NumEle].zmin))
         {
            DummyY[i+MD->NumEle]=MD->Ele[i-MD->NumEle].zmax - MD->Ele[i-MD->NumEle].zmin;
        }
        if(DummyY[i+MD->NumEle]+DummyY[i]>(MD->Ele[i-MD->NumEle].zmax-MD->Ele[i-MD->NumEle].zmin))
        {
            if((DummyY[i]<(MD->Ele[i-MD->NumEle].zmax-MD->Ele[i-MD->NumEle].zmin)))
            {
                DummyY[i]=1.0*((MD->Ele[i-MD->NumEle].zmax-MD->Ele[i-MD->NumEle].zmin)-DummyY[i+MD->NumEle]);
            }
            else //Only possible if DummyY UnSat = Zmax-Zmin : Check to see if this really happens in natural conditions?
            {
                //DummyY[i]=MD->Ele[i-MD->NumEle].zmax - MD->Ele[i-MD->NumEle].zmin;
                DummyY[i+MD->NumEle]=0.0;
            }
        }
      }

    /* Lateral Flux Calculation between Triangular elements Follows  */
    for(i=0; i<MD->NumEle; i++){
         for(j=0; j<3; j++){
              if(MD->Ele[i].nabr[j] > 0){
              /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
                   if(MD->Ele[MD->Ele[i].nabr[j]-1].zmin>MD->Ele[i].zmin){
                        if(MD->Ele[MD->Ele[i].nabr[j]-1].zmin>MD->Ele[i].zmin+DummyY[i+2*MD->NumEle]){
                             Avg_Y_Sub=DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]/2;
                        }
                        else{
                             Avg_Y_Sub=(DummyY[i+2*MD->NumEle]+MD->Ele[i].zmin-MD->Ele[MD->Ele[i].nabr[j]-1].zmin+DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle])/2;
                        }
                   }
                   else{
                        if(MD->Ele[i].zmin>MD->Ele[MD->Ele[i].nabr[j]-1].zmin+DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]){
                             Avg_Y_Sub=DummyY[i+2*MD->NumEle]/2;
                        }
                        else{
                             Avg_Y_Sub=(DummyY[i+2*MD->NumEle]+DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]+MD->Ele[MD->Ele[i].nabr[j]-1].zmin-MD->Ele[i].zmin)/2;
                        }
                   }
                   Dif_Y_Sub = (DummyY[i+2*MD->NumEle] + MD->Ele[i].zmin) - (DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle] + MD->Ele[MD->Ele[i].nabr[j]-1].zmin);
                   Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x), 2) + pow((MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y), 2));
                   Avg_Ksat = (MD->Ele[i].Ksat + MD->Ele[MD->Ele[i].nabr[j] - 1].Ksat)/2.0;
                   Grad_Y_Sub = Dif_Y_Sub/Distance;
                      /* #? see the source of -1 in the nabr[j]-1 */
                   /* take care of macropore effect */
                   if((MD->Ele[i].zmax-MD->Ele[i].zmin)-DummyY[i+2*MD->NumEle]-DummyY[i+MD->NumEle]<=0){
                        elemSatn=1.0;
                   }
                   else{
                        elemSatn = DummyY[i+2*MD->NumEle]/(MD->Ele[i].zmax-MD->Ele[i].zmin);   /*  Will have to change this for other formulation */
                   }
                   if (MD->Soil[(MD->Ele[i].soil-1)].Macropore == 0){
                        mp_factor = 1;
                   }
                   else if (MD->Soil[(MD->Ele[i].soil-1)].Macropore == 1){
                        if((elemSatn>=sat_THRESH)&&(DummyY[i]>ovl_THRESH_H)){
                             temp1=1.0+mp_MULTFH*(elemSatn-sat_THRESH)/(1-sat_THRESH);
                             temp1=temp1*mpArea_CALIB+1*(1-mpArea_CALIB);
                        }
                        else{
                             temp1 = 1.0;
                        }

                        if((MD->Ele[MD->Ele[i].nabr[j] - 1].zmax-MD->Ele[MD->Ele[i].nabr[j] - 1].zmin)-DummyY[MD->Ele[i].nabr[j] - 1+2*MD->NumEle]-DummyY[MD->Ele[i].nabr[j] - 1+MD->NumEle]<=0){
                             elemSatn=1.0;
                        }
                        else{
                             elemSatn = DummyY[MD->Ele[i].nabr[j] - 1+2*MD->NumEle]/(MD->Ele[MD->Ele[i].nabr[j] - 1].zmax-MD->Ele[MD->Ele[i].nabr[j] - 1].zmin);   /*  Will have to change this for other formulation */
                        }
                        if((elemSatn>=sat_THRESH)&&(DummyY[MD->Ele[i].nabr[j] - 1]>ovl_THRESH_H)){
                             temp2=1.0+mp_MULTFH*(elemSatn-sat_THRESH)/(1-sat_THRESH);
                             temp2=temp2*mpArea_CALIB+1*(1-mpArea_CALIB);
                        }
                        else{
                             temp2 = 1.0;
                        }

                        mp_factor = (temp1 + temp2)/2.0;
                   }

                   /* groundwater flow modeled by Darcy's law */
                   MD->FluxSub[i][j] = mp_factor*Kh_CALIB*Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*MD->Ele[i].edge[j];

                   /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */

                   if(DummyY[i + 2*MD->NumEle] <= 0 && MD->FluxSub[i][j] > 0){
                        MD->FluxSub[i][j] = 0;
                   }
                   if(DummyY[MD->Ele[i].nabr[j] - 1 + 2*MD->NumEle] <= 0 && MD->FluxSub[i][j] < 0){
                        MD->FluxSub[i][j] = 0;
                   }

                   /* Saturation check */

                   if((DummyY[i + 2*MD->NumEle] >= (MD->Ele[i].zmax - MD->Ele[i].zmin)) && MD->FluxSub[i][j] < 0){
                        MD->FluxSub[i][j] = 0;
                   }
                   if((DummyY[MD->Ele[i].nabr[j] - 1 + 2*MD->NumEle] >= (MD->Ele[MD->Ele[i].nabr[j] - 1].zmax - MD->Ele[MD->Ele[i].nabr[j] - 1].zmin)) && MD->FluxSub[i][j] >0){
                        MD->FluxSub[i][j] = 0;
                   }

                   /* Surface Lateral Flux Calculation between Triangular elements Follows */
                   if(MD->Ele[MD->Ele[i].nabr[j] - 1].zmax>MD->Ele[i].zmax){
                        if(MD->Ele[MD->Ele[i].nabr[j] - 1].zmax>MD->Ele[i].zmax+DummyY[i]){
                             Avg_Y_Surf=DummyY[MD->Ele[i].nabr[j]-1]/2;
                        }
                        else{
                             Avg_Y_Surf=(DummyY[i]+MD->Ele[i].zmax-MD->Ele[MD->Ele[i].nabr[j] - 1].zmax+DummyY[MD->Ele[i].nabr[j]-1])/2;
                        }
                   }
                   else{
                        if(MD->Ele[i].zmax>MD->Ele[MD->Ele[i].nabr[j] - 1].zmax+DummyY[MD->Ele[i].nabr[j]-1]){
                             Avg_Y_Surf=DummyY[i]/2;
                        }
                        else{
                             Avg_Y_Surf=(DummyY[i]+DummyY[MD->Ele[i].nabr[j]-1]+MD->Ele[MD->Ele[i].nabr[j] - 1].zmax-MD->Ele[i].zmax)/2;
                        }
                   }
                   Dif_Y_Surf = (DummyY[i] + MD->Ele[i].zmax) - (DummyY[MD->Ele[i].nabr[j] - 1] + MD->Ele[MD->Ele[i].nabr[j] - 1].zmax);
                   Grad_Y_Surf = Dif_Y_Surf/Distance;
                   Avg_Sf = (MD->Ele[i].Sf + MD->Ele[MD->Ele[i].nabr[j] - 1].Sf)/2.0;
                   Avg_Rough = 0.5*(MD->Ele[i].Rough + MD->Ele[MD->Ele[i].nabr[j] - 1].Rough);
                   CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];

                   OverlandFlow(MD->FluxSurf,i,j,MD->SurfMode, Avg_Y_Surf,Grad_Y_Surf,Avg_Sf,Alfa,Beta,CrossA,Avg_Rough,1,1);

                    if(isnan(MD->FluxSurf[i][j])==1){
                        printf("\n1: %f %d %d %lf %lf %lf",t,MD->Ele[i].index,MD->Ele[MD->Ele[i].nabr[j]-1].index,DummyY[i],DummyY[MD->Ele[i].nabr[j]-1],MD->FluxSurf[i][j]);
                        getchar();
                   }
                   /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */

                   if(DummyY[i] <= 0 && MD->FluxSurf[i][j] > 0){
                        MD->FluxSurf[i][j] = 0;
                   }
                   if(DummyY[MD->Ele[i].nabr[j] - 1] <= 0 && MD->FluxSurf[i][j] < 0){
                        MD->FluxSurf[i][j] = 0;
                   }
              }
              else{
                   /* Handle boundary conditions for elements. No flow (natural) boundary condition is default */
                   if(MD->Ele[i].BC == 0){
                        MD->FluxSurf[i][j] = 0;
                        MD->FluxSub[i][j] = 0;
                   }
                   else{
                        //MD->FluxSurf[i][j] = 0;

                        /* this part of code has not been tested! */
                        if(MD->Ele[i].BC > 0){         // Dirichlet BC
                loc_bcEle = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC)-1], t);

                                if( loc_bcEle < MD->Ele[i].zmin ){
                                        Avg_Y_Surf = (DummyY[i])/2;
                                        Avg_Y_Sub  = (DummyY[i+2*MD->NumEle])/2;
                                }
                                else if(loc_bcEle < MD->Ele[i].zmax){
                                        Avg_Y_Surf = (DummyY[i])/2;
                                        Avg_Y_Sub  = (loc_bcEle - MD->Ele[i].zmin + DummyY[i+2*MD->NumEle])/2.0;
                                }
                                else{
                                        Avg_Y_Surf = (DummyY[i] + loc_bcEle - MD->Ele[i].zmax)/2;
                                        Avg_Y_Sub  = MD->Ele[i].zmax - MD->Ele[i].zmin;
                                }
                                Distance = sqrt(pow(MD->Ele[i].edge[0]*MD->Ele[i].edge[1]*MD->Ele[i].edge[2]/(4*MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j]/2, 2));

                                //Sub
                                Dif_Y_Sub = DummyY[i+2*MD->NumEle] + MD->Ele[i].zmin - loc_bcEle;
                                Avg_Ksat = MD->Ele[i].Ksat;
                                Grad_Y_Sub = Dif_Y_Sub/Distance;
                                MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];

                                //Surf
                                Dif_Y_Surf = DummyY[i] + MD->Ele[i].zmax - loc_bcEle;
                                Grad_Y_Surf = Dif_Y_Surf / Distance;
                                MD->FluxSurf[i][j] = 0.1*(Grad_Y_Surf>0?1:-1)* (Avg_Y_Sub * MD->Ele[i].edge[j] ) * (pow(pow(Avg_Y_Surf, 1.0/3.0),2)/(MD->Ele[i].Rough)) * sqrt((Grad_Y_Surf>0?1:-1)*Grad_Y_Surf);

                                MD->FluxSurf[i][j] = sqrt(DummyY[i]/Distance)*pow(DummyY[i],2.0/3.0)*(DummyY[i]/2)*MD->Ele[i].edge[j]/MD->Ele[i].Rough;
                        }

                                 /* Note the distance calculated here is the distance between circumcenter of the triangle and the edge on which BDD. condition is defined
                             Distance = sqrt(pow(MD->Ele[i].edge[0]*MD->Ele[i].edge[1]*MD->Ele[i].edge[2]/(4*MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j]/2, 2));
                             Avg_Ksat = MD->Ele[i].Ksat;
                             Grad_Y_Sub = Dif_Y_Sub/Distance;

                             MD->FluxSub[i][j] = Avg_Ksat*Grad_Y_Sub*Avg_Y_Sub*MD->Ele[i].edge[j];
                        } */


                        else{                          // Nuemann BC
                             MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC)-1+MD->Num1BC], t);
                        }
            /*
            Distance = sqrt(pow(MD->Ele[i].edge[0]*MD->Ele[i].edge[1]*MD->Ele[i].edge[2]/(4*MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j]/2, 2));
                        MD->FluxSub[i][j]=MD->Ele[i].Ksat*(DummyY[i+2*MD->NumEle]/Distance)*(DummyY[i+2*MD->NumEle]/2)*MD->Ele[i].edge[j];
            */
                        /* zero-depth-gradient boundary conditions */
            /*
                            Dif_Y_Surf=DummyY[i]/Distance;
                        Avg_Rough = MD->Ele[i].Rough;
                        Avg_Y_Riv = DummyY[i]/2;
                        MD->FluxSurf[i][j] = sqrt(Dif_Y_Surf)*pow(DummyY[i],2.0/3.0)*Avg_Y_Riv*MD->Ele[i].edge[j]/Avg_Rough;
            */
                   }

                   /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */
                   if(DummyY[i + 2*MD->NumEle] <= 0 && MD->FluxSub[i][j] > 0){
                        MD->FluxSub[i][j] = 0;
                   }
              }
         }
         /************************************************************************************************/

         /**************************************Evaporation from ground*********************************************************/
         if((MD->Ele[i].zmax-MD->Ele[i].zmin)-DummyY[i+2*MD->NumEle]-DummyY[i+MD->NumEle]<=0){
              elemSatn=1.0;
            }
         else{
              //elemSatn = 0.5*(1-cos(3.14*(DummyY[i+MD->NumEle]/(MD->Ele[i].zmax-MD->Ele[i].zmin-DummyY[i+2*MD->NumEle]))));    /*  Will have to change this for other formulation */
              elemSatn = (DummyY[i+2*MD->NumEle]+DummyY[i+MD->NumEle]>(MD->Ele[i].zmax-MD->Ele[i].zmin)-MD->Ele[i].RzD-rzd_CALIB)? 0.5*(1-cos(3.14*(DummyY[i+2*MD->NumEle]/(MD->Ele[i].zmax-MD->Ele[i].zmin)))):0;
         }
         Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn-1], t);
         //G = Interpolation(&MD->TSD_G[MD->Ele[i].G-1], t);
         T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp-1], t);
         Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel-1], t);
         RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity-1], t);
         VP = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure-1], t);
         P = 101.325*pow(10,3)*pow((293-0.0065*MD->Ele[i].zmax)/293,5.26);
         Delta = 2503*pow(10,3)*exp(17.27*T/(T+237.3))/(pow(237.3 + T, 2));
         Gamma = P*1.0035*0.92/(0.622*2441);
         LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC-1], t);
         zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
         cnpy_h = zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));

         /*if(LAI<2.85)
        {
        rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
        }
        else
        {
        rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h));
        }
        */
        rl=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
        r_a = log(MD->Ele[i].windH/rl)*log(10*MD->Ele[i].windH/rl)/(Vel*0.16);

 /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DELETE 0.01 BELOW $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

        MD->EleET[i][2] = (1-MD->Ele[i].VegFrac)*(Rn*(1-MD->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a))/(1000.0*2441000.0*(Delta+Gamma));
        MD->EleET[i][2]=et2_CALIB*MD->EleET[i][2];

        if(LAI>0.0){
             Rmax = 5000.0/(24*3600);        /* Unit day_per_m */
             f_r= 1.1*Rn*(1-exp(-LAI))/(MD->Ele[i].Rs_ref*LAI);
             alpha_r= (1+f_r)/(1+(MD->Ele[i].Rmin/Rmax));
             eta_s= 1- 0.0016*(pow((24.85-T),2));
             beta_s= elemSatn>EPSILON/1000.0?elemSatn:EPSILON/1000.0;

             r_s=(MD->Ele[i].Rmin*alpha_r/(beta_s*LAI*pow(eta_s,4)))> Rmax?Rmax:(MD->Ele[i].Rmin*alpha_r/(beta_s*LAI*pow(eta_s,4)));

             MD->EleET[i][1] = (LAI/MD->Ele[i].LAImax)*MD->Ele[i].VegFrac*(1-pow((MD->EleIS[i]<0?0:MD->EleIS[i])/MD->EleISmax[i],2.0/3))*(Rn*(1-MD->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a))/(1000*2441000.0*(Delta+Gamma*(1+r_s/r_a)));
        }
        else{
             MD->EleET[i][1] =0.0;
        }


        /**********************************************************************************************************************/

        MD->EleET[i][1]=et1_CALIB*MD->EleET[i][1];


        AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;

        if(DummyY[i+2*MD->NumEle]+DummyY[i+MD->NumEle] >= AquiferDepth){
             Deficit = 0;
             MD->EleVic[i] = 0;
        }
        else if(DummyY[i+2*MD->NumEle]+DummyY[i+MD->NumEle]>AquiferDepth-MD->Ele[i].RzD){
               Deficit = AquiferDepth - DummyY[i+2*MD->NumEle];
               MD->EleVic[i] = MD->Ele[i].Ksat*(1+(DummyY[i]/MD->Ele[i].RzD));
        }
        else{
             Deficit = AquiferDepth - DummyY[i+2*MD->NumEle];
             MD->EleVic[i] = MD->Ele[i].Ksat*(1+(DummyY[i]-log((DummyY[i+MD->NumEle]+EPSILON/10000.0)/(Deficit-MD->Ele[i].RzD+EPSILON/10000.0))/MD->Ele[i].Alpha)/(MD->Ele[i].RzD));/* Interpolation(&MD->TSD_Inc[MD->Soil[(MD->Ele[i].soil-1)].Inf-1], t);*/    /* ## Take care of this using saturation rate instead of */
        }

        MD->EleVic[i]=MD->EleVic[i]/Vic_CALIB;

        if(DummyY[i+MD->NumEle] < Deficit){                                     /*  ## redo this condition with only phi and z */
             if(DummyY[i] > 0){
             /***************************Addn***************************/
                  if((MD->Ele[i].zmax-MD->Ele[i].zmin)-DummyY[i+2*MD->NumEle]-DummyY[i+MD->NumEle]<=0){
                       elemSatn=1.0;
                  }
                  else{
                       elemSatn = DummyY[i+2*MD->NumEle]/(MD->Ele[i].zmax-MD->Ele[i].zmin);   /*  Will have to change this for other formulation */
                  }

                  if (MD->Soil[(MD->Ele[i].soil-1)].Macropore == 0){
                       mp_factor = 1.0;
                  }
                  else if (MD->Soil[(MD->Ele[i].soil-1)].Macropore == 1){
                       if((elemSatn>=sat_THRESH)&&(DummyY[i]>ovl_THRESH_V)){
                            mp_factor=1.0+mp_MULTFV*(elemSatn-sat_THRESH)/(1-sat_THRESH);
                            mp_factor=mp_factor*mpArea_CALIB+1*(1-mpArea_CALIB);
                       }
                       else{
                            mp_factor = 1.0;
                       }
                  }

                  /*********************************************************/
                  if((MD->Ele[i].zmax-MD->Ele[i].zmin)-DummyY[i+2*MD->NumEle]-DummyY[i+MD->NumEle]<=0){
                       elemSatn=1.0;
                  }
                  else{
                       elemSatn = (DummyY[i+2*MD->NumEle]+DummyY[i+MD->NumEle]>(MD->Ele[i].zmax-MD->Ele[i].zmin)-MD->Ele[i].RzD-rzd_CALIB)? 0.5*(1-cos(3.14*(DummyY[i+2*MD->NumEle]/(MD->Ele[i].zmax-MD->Ele[i].zmin)))):0;
                  }

                  DummyDY[i] = MD->EleNetPrep[i] - mp_factor*MD->EleVic[i] - MD->EleET[i][2];
                  DummyDY[i+MD->NumEle] = (1.0-mpArea_CALIB)*MD->EleVic[i];                    /* ##Clean DY variables to some dummy for defining right hand side */
                  DummyDY[i+2*MD->NumEle] = DummyDY[i+2*MD->NumEle]+(mp_factor-1*(1-mpArea_CALIB))*MD->EleVic[i];
             }
             else if((MD->EleNetPrep[i] > MD->EleVic[i]+MD->EleET[i][2])){
                  DummyDY[i] = MD->EleNetPrep[i] - (MD->EleVic[i]+MD->EleET[i][2]);
                  DummyDY[i+MD->NumEle] = MD->EleVic[i];
             }
             else if((MD->EleNetPrep[i] < MD->EleVic[i]+MD->EleET[i][2])){       /* ##Think abt putting <= here */

                  if(MD->EleNetPrep[i] < MD->EleVic[i]){
                       DummyDY[i] = 0;
                         DummyDY[i+MD->NumEle] = MD->EleNetPrep[i]-elemSatn*MD->EleET[i][2];
                //BHATT
                MD->EleET[i][2]=elemSatn*MD->EleET[i][2];
                  }
                  else{
                         DummyDY[i] = 0;
                         DummyDY[i+MD->NumEle] = MD->EleVic[i];
                  }
             }
             else{
                  DummyDY[i] = 0;
                  DummyDY[i+MD->NumEle] = MD->EleNetPrep[i]-elemSatn*MD->EleET[i][2];
            //BHATT
            MD->EleET[i][2]=elemSatn*MD->EleET[i][2];
             }
        }
        else{                                    /*  ## redo this condition */

             /* The reason of this happening is not clear, therefore the treatment is not secure
             Fortranately, this will not happen in most cases. */
             if(DummyY[i]>0){
                  DummyDY[i] = MD->EleNetPrep[i]-MD->EleET[i][2];
                   DummyDY[i+2*MD->NumEle] = 0;
             }
              else{
                  if(MD->EleNetPrep[i]>0){
                       if(MD->EleNetPrep[i]>MD->EleET[i][2]){
                            DummyDY[i] =MD->EleNetPrep[i]-MD->EleET[i][2];
                         }
                          else{
                               DummyDY[i] = 0;
                          }
                          DummyDY[i+2*MD->NumEle] = 0;
                  }
                   else{
                        DummyDY[i] = 0;
                        DummyDY[i+2*MD->NumEle] = 0-elemSatn*MD->EleET[i][2];
                    //BHATT
                    MD->EleET[i][2]=elemSatn*MD->EleET[i][2];
                  }
             }
        }

        if(DummyY[i+2*MD->NumEle]>=AquiferDepth-MD->Ele[i].RzD){
             DummyDY[i+2*MD->NumEle]=DummyDY[i+2*MD->NumEle] - MD->EleET[i][1];
        }
        else{
              DummyDY[i+MD->NumEle]=DummyDY[i+MD->NumEle] - MD->EleET[i][1];
        }

        /*for(j=0; j<3; j++)
        {
        DummyDY[i] =  DummyDY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
        }

        for(j=0; j<3; j++)
        {
         DummyDY[i+2*MD->NumEle] = DummyDY[i+2*MD->NumEle] - MD->FluxSub[i][j]/MD->Ele[i].area;
        } */
         /************************************************************************************************/
    }


    /* initialize river flux */
    for(i=0; i<MD->NumRiv; i++){
         for(j=0; j<6; j++){
              MD->FluxRiv[i][j] = 0;
         }
    }

    /* Lateral Flux Calculation between River-River and River-Triangular elements Follows */
    for(i=0; i<MD->NumRiv; i++){

         /* Note: the ordering of river segment in input file has to be: from UP to DOWN */
         TotalY_Riv = DummyY[i + 3*MD->NumEle] + MD->Riv[i].zmin;
         Perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[i + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
         /*    if(DummyY[10 + 3*MD->NumEle]>0)
         {
               printf("\n%lf %e %e %e %e area",t,DummyY[10 + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[10].shape - 1].coeff,Perem,CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[10].shape - 1].interpOrd,DummyY[10 + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[10].shape - 1].coeff,1));
               getchar();
         }
         */
         /* Lateral Flux Calculation between River-River element Follows */
         if(MD->Riv[i].down > 0){
              TotalY_Riv_down = DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle] + MD->Riv[MD->Riv[i].down - 1].zmin;
              Perem_down = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].coeff,2);
              Avg_Perem = (Perem + Perem_down)/2.0;    /* Avg perimeter */
              if(MD->Riv[MD->Riv[i].down - 1].zmin>MD->Riv[i].zmin){
                   if(MD->Riv[MD->Riv[i].down - 1].zmin>MD->Riv[i].zmin+DummyY[i + 3*MD->NumEle]){
                        Avg_Y_Riv=DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle]/2;
                   }
                   else{
                        Avg_Y_Riv=(DummyY[i + 3*MD->NumEle]+MD->Riv[i].zmin-MD->Riv[MD->Riv[i].down - 1].zmin+DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle])/2;
                   }
              }
              else{
                   if(MD->Riv[i].zmin>MD->Riv[MD->Riv[i].down - 1].zmin+DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle]){
                        Avg_Y_Riv=DummyY[i + 3*MD->NumEle]/2;
                   }
                   else{
                        Avg_Y_Riv=(DummyY[i + 3*MD->NumEle]+DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle]+MD->Riv[MD->Riv[i].down - 1].zmin-MD->Riv[i].zmin)/2;
                   }
              }
             Avg_Rough = (MD->Riv_Mat[MD->Riv[i].material - 1].Rough + MD->Riv_Mat[MD->Riv[MD->Riv[i].down - 1].material-1].Rough)/2.0;
             //Distance = sqrt(pow(MD->Riv[i].x - MD->Riv[MD->Riv[i].down - 1].x, 2) + pow(MD->Riv[i].y - MD->Riv[MD->Riv[i].down - 1].y, 2));
             Distance = (MD->Riv[i].Length+MD->Riv[MD->Riv[i].down - 1].Length)/2;

             Dif_Y_Riv = (TotalY_Riv - TotalY_Riv_down)/Distance;
             Avg_Sf = (MD->Riv_Mat[MD->Riv[i].material - 1].Sf + MD->Riv_Mat[MD->Riv[MD->Riv[i].down - 1].material-1].Sf)/2.0;
             /*CrossA = 0.5*(CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[i + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,1)+CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd,DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].coeff,1));
             */
             CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,Avg_Y_Riv,MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,1);
             OverlandFlow(MD->FluxRiv,i,1,MD->RivMode, Avg_Y_Riv,Dif_Y_Riv,Avg_Sf,Alfa,Beta,CrossA,Avg_Rough,0,Avg_Perem);

             /* Correction is being done in flux terms which can be > 0 even when there is no source water level present */
             if(DummyY[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] > 0){
                  MD->FluxRiv[i][1] = 0;
             }
             else if(DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] < 0){
                  MD->FluxRiv[i][1] = 0;
             }

             /* accumulate to get in-flow for down segments: [0] for inflow, [1] for outflow */
             MD->FluxRiv[MD->Riv[i].down - 1][0] = MD->FluxRiv[MD->Riv[i].down - 1][0] + MD->FluxRiv[i][1];
         }
         /* Here the if statement ends for if there is a downstream river segment or not */
         else{
              switch(MD->Riv[i].down){
                   case -1:

                        /* Dirichlet boundary condition */
                        TotalY_Riv_down = Interpolation(&MD->TSD_Riv[(MD->Riv[i].BC)-1], t) + MD->Node[MD->Riv[i].ToNode-1].zmin + MD->Riv_Shape[MD->Riv[i].shape-1].bed;
                        Distance = (MD->Riv[i].Length)*0.5;
                        Dif_Y_Riv = (TotalY_Riv - TotalY_Riv_down)/Distance;
                        Avg_Sf = MD->Riv_Mat[MD->Riv[i].material - 1].Sf;
                        Avg_Rough = MD->Riv_Mat[MD->Riv[i].material-1].Rough;
                        Avg_Y_Riv = DummyY[i + 3*MD->NumEle];
                        Avg_Perem = Perem;
                        CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[i + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,1);

                        OverlandFlow(MD->FluxRiv,i,1,MD->RivMode, Avg_Y_Riv,Dif_Y_Riv,Avg_Sf,Alfa,Beta,CrossA,Avg_Rough,0,Avg_Perem);

                        break;

                   case -2:

                        /* Neumann boundary condition */
                        MD->FluxRiv[i][1] = Interpolation(&MD->TSD_Riv[MD->Riv[i].BC-1], t);
                        break;

                   case -3:

                        /* zero-depth-gradient boundary conditions */
                        Distance = (MD->Riv[i].Length)*0.5;
                        //Dif_Y_Riv = (MD->Riv[i].zmin - (MD->Node[MD->Riv[i].ToNode-1].zmax -MD->Riv[i].depth+ MD->Riv_Shape[MD->Riv[i].shape-1].bed))/Distance;
                          Dif_Y_Riv=0.1/Distance;
                        Avg_Rough = MD->Riv_Mat[MD->Riv[i].material-1].Rough;
                        Avg_Y_Riv = DummyY[i + 3*MD->NumEle];
                        Avg_Perem = Perem;
                        CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[i + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,1);
                        MD->FluxRiv[i][1] = sqrt(Dif_Y_Riv)*CrossA*(Perem>0?pow(CrossA/Perem,2.0/3.0):0)/Avg_Rough;
                        break;
                        /* #? How is critical dept being defined */
                   case -4:

                        /* Critical Depth boundary conditions */
                        CrossA = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[i + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,1);
                        MD->FluxRiv[i][1] = CrossA*sqrt(GRAV*DummyY[i + 3*MD->NumEle]);
                        break;

                   default:
                        printf("Fatal Error: River Routing Boundary Condition Type Is Wrong!");
                        exit(1);
              }

              /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */
              if(DummyY[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][1] > 0){
                   MD->FluxRiv[i][1] = 0;
              }
              /*      else if(DummyY[MD->Riv[i].down - 1 + 2*MD->NumEle] <= 0 && MD->FluxRiv[i][1] < 0)
              {
                      MD->FluxRiv[i][1] = 0;
                       }  */            /* ######## Does this need to be commented */

              /* #? Have to take care of this before it can be implemented on ny basin other than ShaleHIlls */
              /* need some work to take care of multiple output Q */
              MD->Q = MD->FluxRiv[i][1];

         }

         /* Lateral Surface Flux Calculation between River-Triangular element Follows */
         if (MD->Riv[i].LeftEle > 0){
              OLflowFromEleToRiv(DummyY[MD->Riv[i].LeftEle - 1],MD->Ele[MD->Riv[i].LeftEle - 1].zmax,MD->Riv[i].x,MD->Ele[MD->Riv[i].LeftEle - 1].x,MD->Riv[i].y,MD->Ele[MD->Riv[i].LeftEle - 1].y,MD->Riv_Mat[MD->Riv[i].material-1].Cwr, MD->Riv[i].zmax,TotalY_Riv,MD->FluxRiv,i,2,MD->Riv[i].Length);

              /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */
              if(DummyY[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][2] > 0){
                   MD->FluxRiv[i][2] = 0;
              }

              if(DummyY[MD->Riv[i].LeftEle - 1] <= 0 && MD->FluxRiv[i][2] < 0){
                   MD->FluxRiv[i][2] = 0;
              }

              /* replace overland flux item */
              for(j=0; j < 3; j++){
                   /* this may cause trouble when there are 2 boundary side and one is not neutral.*/
                   if(MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j] == MD->Riv[i].RightEle){
                        MD->FluxSurf[MD->Riv[i].LeftEle - 1][j] = -MD->FluxRiv[i][2];
                   }
              }
         }

         if (MD->Riv[i].RightEle > 0){
              OLflowFromEleToRiv(DummyY[MD->Riv[i].RightEle - 1],MD->Ele[MD->Riv[i].RightEle - 1].zmax,MD->Riv[i].x,MD->Ele[MD->Riv[i].RightEle - 1].x,MD->Riv[i].y,MD->Ele[MD->Riv[i].RightEle - 1].y,MD->Riv_Mat[MD->Riv[i].material-1].Cwr, MD->Riv[i].zmax,TotalY_Riv,MD->FluxRiv,i,3,MD->Riv[i].Length);

              /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */
              if(DummyY[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][3] > 0){
                   MD->FluxRiv[i][3] = 0;
              }

              if(DummyY[MD->Riv[i].RightEle - 1] <= 0 && MD->FluxRiv[i][3] < 0){
                   MD->FluxRiv[i][3] = 0;
              }

              /* replace overland flux item */
              for(j=0; j < 3; j++){
                   if(MD->Ele[MD->Riv[i].RightEle - 1].nabr[j] == MD->Riv[i].LeftEle){
                        MD->FluxSurf[MD->Riv[i].RightEle - 1][j] = -MD->FluxRiv[i][3];
                        break;
                   }
              }
         }
         /* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
         if (MD->Riv[i].LeftEle > 0){
              /* Left Neighbor Groundwater Head */
              /*********************TRANSFER THIS INSIDE FUNC*****************************/
              if((MD->Ele[MD->Riv[i].LeftEle - 1].zmax-MD->Ele[MD->Riv[i].LeftEle - 1].zmin)-DummyY[MD->Riv[i].LeftEle - 1+2*MD->NumEle]-DummyY[MD->Riv[i].LeftEle - 1+MD->NumEle]<=0){
                   elemSatn=1.0;
              }
              else{
                   elemSatn = DummyY[MD->Riv[i].LeftEle - 1+2*MD->NumEle]/(MD->Ele[MD->Riv[i].LeftEle - 1].zmax-MD->Ele[MD->Riv[i].LeftEle - 1].zmin);   /*  Will have to change this for other formulation */
              }

              if (MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].Macropore == 0){
                   mp_factor = 1;
              }
              else if (MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].Macropore == 1){
                   if((elemSatn>=sat_THRESH)&&(DummyY[MD->Riv[i].LeftEle - 1]>ovl_THRESH_H)){
                        mp_factor=1.0+mp_MULTFH*(elemSatn-sat_THRESH)/(1-sat_THRESH);
                        mp_factor=mp_factor*mpArea_CALIB+1*(1-mpArea_CALIB);
                   }
                   else{
                        mp_factor = 1.0;
                   }
              }
              /**************************************************************************/


              if( (MD->Riv[i].zmin + DummyY[i+3*MD->NumEle] - (DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle - 1].zmin))>0 ){
                      loc_perem=Perem;
              }
              else{
                    if(MD->Ele[MD->Riv[i].LeftEle - 1].zmin < MD->Riv[i].zmin){
                          if( (DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle - 1].zmin - MD->Riv[i].zmin) > MD->Riv[i].depth ){
                              loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);

                        }
                          else{
                              loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,(DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle - 1].zmin - MD->Riv[i].zmin),MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
                        }
                    }
                    else{
                        if( DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] > MD->Riv[i].depth){
                            loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
                        }
                        else{
                            loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
                        }
                    }
              }

              GWflowFromEleToRiv(DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle],MD->Ele[MD->Riv[i].LeftEle - 1].zmax,MD->Ele[MD->Riv[i].LeftEle - 1].zmin,MD->Riv[i].x,MD->Ele[MD->Riv[i].LeftEle - 1].x,MD->Riv[i].y,MD->Ele[MD->Riv[i].LeftEle - 1].y,MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].Macropore,DummyY[i+3*MD->NumEle],TotalY_Riv,MD->FluxRiv,i,4,MD->Riv[i].Length,MD->Soil[(MD->Ele[MD->Riv[i].LeftEle - 1].soil-1)].base,mp_factor,loc_perem,MD->Ele[MD->Riv[i].LeftEle - 1].Ksat,MD->Ele[MD->Riv[i].LeftEle-1].RzD); /* delete replace Wid by 0.5*avg_perim */

              /* Saturation check */
              if((DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] >= MD->Ele[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle].zmax-MD->Ele[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle].zmin) && MD->FluxRiv[i][4] > 0){
                   MD->FluxRiv[i][4]  = 0;
              }


              /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */
              if(DummyY[i+3*MD->NumEle] <= 0 && MD->FluxRiv[i][4] > 0){
                   MD->FluxRiv[i][4] = 0;
              }
              if(DummyY[MD->Riv[i].LeftEle - 1 + 2*MD->NumEle] <= 0 && MD->FluxRiv[i][4] < 0){
                   MD->FluxRiv[i][4] = 0;
              }

              /* modify groundwater flux item */
              for(j=0; j < 3; j++){
                   if(MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j] == MD->Riv[i].RightEle){
                        if(MD->Ele[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1].zmin>MD->Ele[MD->Riv[i].LeftEle - 1].zmin){
                                if(MD->Ele[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1].zmin>MD->Ele[MD->Riv[i].LeftEle - 1].zmin+DummyY[MD->Riv[i].LeftEle - 1+2*MD->NumEle]){
                                   Avg_Y_Sub=DummyY[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1 + 2*MD->NumEle]/2;
                             }
                                else{
                                      Avg_Y_Sub=(DummyY[MD->Riv[i].LeftEle - 1+2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle - 1].zmin-MD->Ele[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1].zmin+DummyY[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1 + 2*MD->NumEle])/2;
                             }
                        }
                           else{
                                if(MD->Ele[MD->Riv[i].LeftEle - 1].zmin>MD->Ele[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1].zmin+DummyY[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1 + 2*MD->NumEle]){
                                   Avg_Y_Sub=DummyY[MD->Riv[i].LeftEle - 1+2*MD->NumEle]/2;
                                }
                                else{
                                      Avg_Y_Sub=(DummyY[MD->Riv[i].LeftEle - 1+2*MD->NumEle]+DummyY[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1 + 2*MD->NumEle]+MD->Ele[MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j]-1].zmin-MD->Ele[MD->Riv[i].LeftEle - 1].zmin)/2;
                                }
                           }
                        Avg_BedDepth = (MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin + MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)*0.5;
                        if(Avg_BedDepth - MD->Riv[i].depth <= 0)
                        {
                            MD->FluxSub[MD->Riv[i].LeftEle - 1][j] = 0.0;
                        }
                        else{
                            if (Avg_Y_Sub > Avg_BedDepth - MD->Riv[i].depth){
                                MD->FluxSub[MD->Riv[i].LeftEle - 1][j] = MD->FluxSub[MD->Riv[i].LeftEle - 1][j] * (Avg_BedDepth - MD->Riv[i].depth) / Avg_Y_Sub;
                            }
                        }
                        break;
                   }
              }
         }

         if (MD->Riv[i].RightEle > 0){
         /************************TRANSFER THIS INSIDE FUNC*******************/
              if((MD->Ele[MD->Riv[i].RightEle - 1].zmax-MD->Ele[MD->Riv[i].RightEle - 1].zmin)-DummyY[MD->Riv[i].RightEle - 1+2*MD->NumEle]-DummyY[MD->Riv[i].RightEle - 1+MD->NumEle]<=0){
                   elemSatn=1.0;
              }
              else{
                   elemSatn = DummyY[MD->Riv[i].RightEle - 1+2*MD->NumEle]/(MD->Ele[MD->Riv[i].RightEle - 1].zmax-MD->Ele[MD->Riv[i].RightEle - 1].zmin);   /*  Will have to change this for other formulation */
              }

              if (MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].Macropore == 0){
                   mp_factor = 1;
              }
              else if (MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].Macropore == 1){
                   if((elemSatn>=sat_THRESH)&&(DummyY[MD->Riv[i].RightEle - 1]>ovl_THRESH_H)){
                        mp_factor=1.0+mp_MULTFH*(elemSatn-sat_THRESH)/(1-sat_THRESH);
                        mp_factor=mp_factor*mpArea_CALIB+1*(1-mpArea_CALIB);
                   }
                   else{
                        mp_factor = 1.0;
                   }
              }
              /*******************************************************************/
              if( (MD->Riv[i].zmin + DummyY[i+3*MD->NumEle] - (DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle - 1].zmin))>0 ){
                      loc_perem=Perem;
              }
              else{
                    if(MD->Ele[MD->Riv[i].RightEle - 1].zmin < MD->Riv[i].zmin){
                          if( (DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].RightEle - 1].zmin - MD->Riv[i].zmin) > MD->Riv[i].depth ){
                              loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);

                        }
                          else{
                              loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,(DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] + MD->Ele[MD->Riv[i].RightEle - 1].zmin - MD->Riv[i].zmin),MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
                        }
                    }
                    else{
                        if( DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] > MD->Riv[i].depth){
                            loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
                        }
                        else{
                            loc_perem = CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,2);
                        }
                    }
              }

              GWflowFromEleToRiv(DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle],MD->Ele[MD->Riv[i].RightEle - 1].zmax,MD->Ele[MD->Riv[i].RightEle - 1].zmin,MD->Riv[i].x,MD->Ele[MD->Riv[i].RightEle - 1].x,MD->Riv[i].y,MD->Ele[MD->Riv[i].RightEle - 1].y,MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].Macropore,DummyY[i+3*MD->NumEle],TotalY_Riv,MD->FluxRiv,i,5,MD->Riv[i].Length,MD->Soil[(MD->Ele[MD->Riv[i].RightEle - 1].soil-1)].base,mp_factor,loc_perem,MD->Ele[MD->Riv[i].RightEle - 1].Ksat,MD->Ele[MD->Riv[i].RightEle-1].RzD); /* delete replace Wid by 0.5*avg_perim */

              /* Saturation check */
              if((DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] >= MD->Ele[MD->Riv[i].RightEle - 1 + 2*MD->NumEle].zmax-MD->Ele[MD->Riv[i].RightEle - 1 + 2*MD->NumEle].zmin) && MD->FluxRiv[i][5] > 0){
                   MD->FluxRiv[i][5]  = 0;
              }

              /*     Correction is being done in flux terms which can be > 0 even when there is no source water level present */
              if(DummyY[i + 3*MD->NumEle] <= 0 && MD->FluxRiv[i][5] > 0){
                   MD->FluxRiv[i][5] = 0;
              }
              if(DummyY[MD->Riv[i].RightEle - 1 + 2*MD->NumEle] <= 0 && MD->FluxRiv[i][5] < 0){
                   MD->FluxRiv[i][5] = 0;
              }

              /* modify groundwater flux item */
              for(j=0; j < 3; j++){
                  if(MD->Ele[MD->Riv[i].RightEle - 1].nabr[j] == MD->Riv[i].LeftEle){
                        if(MD->Ele[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1].zmin>MD->Ele[MD->Riv[i].RightEle - 1].zmin){
                                if(MD->Ele[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1].zmin>MD->Ele[MD->Riv[i].RightEle - 1].zmin+DummyY[MD->Riv[i].RightEle - 1+2*MD->NumEle]){
                                   Avg_Y_Sub=DummyY[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1 + 2*MD->NumEle]/2;
                             }
                                else{
                                      Avg_Y_Sub=(DummyY[MD->Riv[i].RightEle - 1+2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle - 1].zmin-MD->Ele[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1].zmin+DummyY[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1 + 2*MD->NumEle])/2;
                             }
                        }
                           else{
                                if(MD->Ele[MD->Riv[i].RightEle - 1].zmin>MD->Ele[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1].zmin+DummyY[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1 + 2*MD->NumEle]){
                                   Avg_Y_Sub=DummyY[MD->Riv[i].RightEle - 1+2*MD->NumEle]/2;
                                }
                                else{
                                      Avg_Y_Sub=(DummyY[MD->Riv[i].RightEle - 1+2*MD->NumEle]+DummyY[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1 + 2*MD->NumEle]+MD->Ele[MD->Ele[MD->Riv[i].RightEle - 1].nabr[j]-1].zmin-MD->Ele[MD->Riv[i].RightEle - 1].zmin)/2;
                                }
                           }
                          Avg_BedDepth = (MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin + MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)*0.5;
                        if(Avg_BedDepth - MD->Riv[i].depth <= 0){
                            MD->FluxSub[MD->Riv[i].RightEle - 1][j] = 0.0;
                        }
                        else{
                              if (Avg_Y_Sub > Avg_BedDepth - MD->Riv[i].depth){
                                  MD->FluxSub[MD->Riv[i].RightEle - 1][j] = MD->FluxSub[MD->Riv[i].RightEle - 1][j] * (Avg_BedDepth - MD->Riv[i].depth) / Avg_Y_Sub;
                            }
                          }
                        break;
                  }
              }
         }
         /**************************************************************************************************/
         if (MD->Riv[i].LeftEle > 0){
              DummyDY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] = DummyDY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle] + MD->FluxRiv[i][4]/MD->Ele[MD->Riv[i].LeftEle-1].area;
         }
         if (MD->Riv[i].RightEle > 0){
              DummyDY[MD->Riv[i].RightEle-1 + 2*MD->NumEle] = DummyDY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]
                                 + MD->FluxRiv[i][5]/MD->Ele[MD->Riv[i].RightEle-1].area;
         }

         /**************************************************************************************************/
          if((isnan(MD->FluxRiv[i][0])==1)||(isnan(MD->FluxRiv[i][1])==1)||(isnan(MD->FluxRiv[i][2])==1)||(isnan(MD->FluxRiv[i][3])==1)||(isnan(MD->FluxRiv[i][4])==1)||(isnan(MD->FluxRiv[i][5])==1)){
              printf("\n2: %d %d %d %d :%lf %lf: %lf %lf: %lf %lf: %lf %lf %lf %lf %lf %lf",MD->Riv[i].index,MD->Riv[MD->Riv[i].down-1].index,MD->Riv[i].LeftEle,MD->Riv[i].RightEle,DummyY[i+3*MD->NumEle],DummyY[MD->Riv[i].down-1+3*MD->NumEle],DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle],DummyY[MD->Riv[i].LeftEle-1 + 0*MD->NumEle],DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle],DummyY[MD->Riv[i].RightEle-1 + 0*MD->NumEle],MD->FluxRiv[i][0],MD->FluxRiv[i][1],MD->FluxRiv[i][2],MD->FluxRiv[i][3],MD->FluxRiv[i][4],MD->FluxRiv[i][5]);
              getchar();
         }
    }


    for(i=0; i<MD->NumEle; i++){
         for(j=0; j<3; j++){
              DummyDY[i] =  DummyDY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
              DummyDY[i+2*MD->NumEle] = DummyDY[i+2*MD->NumEle] - MD->FluxSub[i][j]/MD->Ele[i].area;
         }
         /*
         if(i==144){
              //printf("\nIn Flux %f %f %f %f %f %f %f",t,DummyY[i],DummyY[i+MD->NumEle],DummyY[i+2*MD->NumEle],DummyDY[i],DummyDY[i+MD->NumEle],DummyDY[i+2*MD->NumEle]);
         }
         */

        //      AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
        //      Deficit = AquiferDepth - DummyY[i+MD->NumEle];
        //      G = MD->Ele[i].Porosity*(1-pow( 1+pow(MD->Ele[i].Alpha*(Deficit+EPSILON), MD->Ele[i].Beta),-(MD->Ele[i].Beta+1)/MD->Ele[i].Beta));


        /* for test case 1 only */
        /*
        G = MD->Ele[i].Porosity;
        GI = 0;
        */
        /*        if(i==160)
        {
                  printf("\nt= %f y =%.9f sy= %.9f",t,DummyY[i],DummyY[i+MD->NumEle]);
        }
        */
        //DummyDY[i+MD->NumEle] = DummyDY[i+MD->NumEle]/MD->Ele[i].Porosity;        /* ## Explicitly define porosity */

        AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
        Deficit = AquiferDepth - DummyY[i+2*MD->NumEle];
        //      PH = 1 - exp(-MD->Ele[i].Alpha*Deficit);

        //     elemSatn = Deficit-DummyY[i+MD->NumEle]>0?DummyY[i+MD->NumEle]/Deficit:(Deficit==0?1.0:0.0);
        if((MD->Ele[i].zmax-MD->Ele[i].zmin)-DummyY[i+2*MD->NumEle]-DummyY[i+MD->NumEle]<=0){
                elemSatn=1.0;
        }
        else{
             elemSatn = 0.5*(1-cos(3.14*(DummyY[i+MD->NumEle]/(MD->Ele[i].zmax-MD->Ele[i].zmin-DummyY[i+2*MD->NumEle]))));    /*  Will have to change this for other formulation */
        }

        MD->Recharge[i] = elemSatn==0.0?0:Rec_CALIB*(-MD->Ele[i].Ksat*elemSatn*AquiferDepth*(1-(2*log(elemSatn)/(AquiferDepth*MD->Ele[i].Alpha)))/((AquiferDepth-DummyY[i+2*MD->NumEle])+DummyY[i+2*MD->NumEle]*elemSatn));

        if(DummyY[i+MD->NumEle]<=0 && MD->Recharge[i] < 0){
             MD->Recharge[i] = 0;
        }
        if(DummyY[i+2*MD->NumEle]<=0 && MD->Recharge[i] > 0){
             MD->Recharge[i] = 0;
        }

        DummyDY[i+MD->NumEle] = DummyDY[i+MD->NumEle] + MD->Recharge[i];
        DummyDY[i+MD->NumEle] = DummyDY[i+MD->NumEle]/MD->Ele[i].Porosity;

        /* Source check */
        /*      if(DummyY[i+MD->NumEle]>=Deficit-EPSILON && DummyDY[i+MD->NumEle]>0)
        {
        DummyDY[i+MD->NumEle] = 0;
        }
        if(DummyY[i+MD->NumEle]<=0 && DummyDY[i+MD->NumEle]<0)
        {
        DummyDY[i+MD->NumEle] = 0;
        }
        */
        DummyDY[i+2*MD->NumEle] = DummyDY[i+2*MD->NumEle]-MD->Recharge[i];
        DummyDY[i+2*MD->NumEle] = DummyDY[i+2*MD->NumEle]/MD->Ele[i].Porosity;

        if(DummyY[i+2*MD->NumEle]>=AquiferDepth && DummyDY[i+2*MD->NumEle]>0){
             DummyDY[i+2*MD->NumEle] = 0;
        }
        if(DummyY[i+2*MD->NumEle]<=0 && DummyDY[i+2*MD->NumEle]<0){
             DummyDY[i+2*MD->NumEle] = 0;
        }
        if(DummyY[i+MD->NumEle]>=AquiferDepth && DummyDY[i+MD->NumEle]>0){  // WHAT?
             DummyDY[i+MD->NumEle] = 0;
        }
        if(DummyY[i+MD->NumEle]<=0 && DummyDY[i+MD->NumEle]<0){
             DummyDY[i+MD->NumEle] = 0;
        }
        if(DummyY[i]<=0 && DummyDY[i]<0){
             DummyDY[i] = 0;
        }
        // if there is a source (well) in it
        /*      if(MD->Ele[i].source > 0)
        {
        DummyDY[i+MD->NumEle] = DummyDY[i+MD->NumEle] - Interpolation(&MD->TSD_Source[MD->Ele[i].source - 1], t)/(MD->Ele[i].Porosity*MD->Ele[i].area);
        }
        */
        // check if it is out of bound

        /*      if(Y[i+MD->NumEle]>AquiferDepth && DummyDY[i+MD->NumEle]>0)
        {
        DummyDY[i+MD->NumEle] = 0;
        }
        if(Y[i+MD->NumEle]<0 && DummyDY[i+MD->NumEle]<0)
        {
        DummyDY[i+MD->NumEle] = 0;
        }
        */
    }
    for(i=0; i<MD->NumRiv; i++){
         DummyDY[i+3*MD->NumEle] = MD->FluxRiv[i][0] - MD->FluxRiv[i][1];
         DummyDY[i+3*MD->NumEle] = DummyDY[i+3*MD->NumEle] - MD->FluxRiv[i][2] - MD->FluxRiv[i][3];
         DummyDY[i+3*MD->NumEle] = DummyDY[i+3*MD->NumEle] - MD->FluxRiv[i][4] - MD->FluxRiv[i][5];
         DummyDY[i+3*MD->NumEle] = DummyDY[i+3*MD->NumEle]/(MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,DummyY[i + 3*MD->NumEle],MD->Riv_Shape[MD->Riv[i].shape - 1].coeff,3)); /* delete derive denominator to be replace by volume */
         if(DummyY[i+3*MD->NumEle]<=0 && DummyDY[i+3*MD->NumEle]<0){
              DummyDY[i+3*MD->NumEle] = 0;
         }
    }
    /*  for(i=0;i<MD->NumEle;i++)
        {
    printf("{");
    for(j=0;j<3;j++)
        {
        printf("%lf ",MD->FluxSub[i][j]);
        }
    printf("},");
    }
    for(i=0;i<MD->NumEle;i++)
    {
    printf("{");
    for(j=0;j<3;j++)
        {
        printf("%lf ",MD->FluxSurf[i][j]);
        }
    printf("},");
    }
    */

    //printf("Flux: %lf %lf %lf %lf %lf %lf\n", MD->FluxSub[5067-2*MD->NumEle][0], MD->FluxSub[5067-2*MD->NumEle][1], MD->FluxSub[5067-2*MD->NumEle][2], MD->Recharge[5067-2*MD->NumEle],MD->FluxRiv[386][4], MD->FluxRiv[386][5]);
    for(i=0; i<3*MD->NumEle+MD->NumRiv; i++){
         /* Y[i]=DummyY[i]; */
           DY[i]=DummyDY[i]/(60.0*24.0);
    //printf("%d\t%lf\n", i, DY[i]);
    }

    free(DummyY);
    free(DummyDY);
    return(0);
}


realtype Interpolation(TSD *Data, realtype t)
//! Function interpolates the data value at time t from a given Time Series
/*! \param Data is the pointer to a time series data
    \param t is the time of simulation
*/
{
    int i, success;
    realtype result;

    //i=0;
    i=Data->iCounter;
    success = 0;
    t=t/(24.0*60.0);
    while(i<Data->length && t>Data->TS[i][0]){
         i++;
    }

    if(i==0){
         /* t is smaller than the 1st node */
         result = Data->TS[i][1];
    }
    else if(i >= Data->length){
         result = Data->TS[i-1][1];
    }
    else{
         result = ((Data->TS[i][0]-t)*Data->TS[i-1][1] + (t-Data->TS[i-1][0])*Data->TS[i][1])/(Data->TS[i][0]-Data->TS[i-1][0]);
         success = 1;
    }

    if(success == 0){
         /*
         printf("\nWarning:  Exterpolation is used ...\n");
         */
    }

    return result;
}


realtype returnVal(realtype rArea, realtype rPerem, realtype eqWid, realtype ap_Bool)
//! Function returns Area, Peremeter or Equivalent Width depending on the ap_Bool identifier
/*! \param rArea is the area
    \param rPerem is the peremeter
    \param eqWid is the equivalent width
    \param ap_Bool is the identifier to return type
*/
{
    if(ap_Bool==1)
    {
        return rArea;
    }
    else if(ap_Bool==2)
    {
        return rPerem;
    }
    else
    {
        return eqWid;
    }
}

/*    Area or Peremeter of River Segment's cross-section    */
realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool)
//! returns Area or Peremeter of a river segment cross-section
/*! \param rivOrder is the interpolation order of the river segment
    \param rivDepth is the depth of water in the river segment
    \param rivCoeff is the interpolation factor of the river segment
    \param a_pBool is identifer for either Area or Peremeter
*/
{
    realtype rivArea, rivPerem, eq_Wid;
    switch(rivOrder)
    {
        case 1:
            rivArea = rivDepth*rivCoeff;
            rivPerem= 2.0*rivDepth+rivCoeff;
            eq_Wid=rivCoeff;
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        case 2:
            rivArea = pow(rivDepth,2)/rivCoeff;
            rivPerem = 2.0*rivDepth*pow(1+pow(rivCoeff,2),0.5)/rivCoeff;
            eq_Wid=2.0*pow(rivDepth+EPSILON,1/(rivOrder-1))/pow(rivCoeff,1/(rivOrder-1));
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        case 3:
            rivArea = 4*pow(rivDepth,1.5)/(3*pow(rivCoeff,0.5));
            rivPerem =(pow(rivDepth*(1+4*rivCoeff*rivDepth)/rivCoeff,0.5))+(log(2*pow(rivCoeff*rivDepth,0.5)+pow(1+4*rivCoeff*rivDepth,0.5))/(2*rivCoeff));
            eq_Wid=2.0*pow(rivDepth+EPSILON,1/(rivOrder-1))/pow(rivCoeff,1/(rivOrder-1));
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        case 4:
            rivArea = 3*pow(rivDepth,4.0/3.0)/(2*pow(rivCoeff,1.0/3.0));
            rivPerem = 2*((pow(rivDepth*(1+9*pow(rivCoeff,2.0/3.0)*rivDepth),0.5)/3)+(log(3*pow(rivCoeff,1.0/3.0)*pow(rivDepth,0.5)+pow(1+9*pow(rivCoeff,2.0/3.0)*rivDepth,0.5))/(9*pow(rivCoeff,1.0/3.0))));
            eq_Wid=2.0*pow(rivDepth+EPSILON,1/(rivOrder-1))/pow(rivCoeff,1/(rivOrder-1));
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        default:
            printf("\n Relevant Values entered are wrong");
            printf("\n Depth: %lf\tCoeff: %lf\tOrder: %d\t");
            return 0;
    }
}

/*    Function to compute surface-flux between elements    */
void OverlandFlow(realtype **flux, int loci, int locj, int surfmode, realtype avg_y, realtype grad_y, realtype avg_sf, realtype alfa, realtype beta, realtype crossA, realtype avg_rough, int eletypeBool, realtype avg_perem) /* delete perimeter should be passed */
//! Computes surface flux across the edge between two elements or two river segments
/*! \param flux is the pointer to FluxSurf or FluxRiv
    \param loci is Element/River Segment Number
    \param locj is the Neighbour Number of the Element / 1 for River outflow
    \param surfmode is identifier to the Surface Flow mode
    \param avg_y is the avarage head between the elements
    \param grad_y is the hydraulic gradient between the elements
    \param avg_sf is the avarage friction slope of the elements
    \param alfa is dummy variable
    \param beta is dummy variable
    \param crossA is average area of cross-section
    \param avg_rough is avarage manning's roughness coefficient
    \param eletypeBool is an identifier to element mode 1: Element 0: River
    \param avg_perem is the avarage wetted perimeter
*/
{
    int locBool;
    float hydRadius;
    /* if surface gradient is not enough to overcome the friction */
    if(fabs(grad_y) <= avg_sf)
    {
         flux[loci][locj] = 0;
    }
    else
    {
         if(grad_y > 0)
         {
              locBool=1;
         }
         else if(grad_y < 0)
         {
              locBool=-1;
         }
         switch(surfmode)
         {
                case 1:
                  if(eletypeBool==1)
                  {
                       /* Kinematic Wave Approximation constitutive relationship: Manning Equation */
                       alfa = sqrt(locBool*grad_y)/avg_rough;
                       beta = pow(avg_y, 2.0/3.0);
                       flux[loci][locj] = locBool*alfa*beta*crossA;
                       break;
                  }
                  else
                  {
                       /*alfa = sqrt(locBool*grad_y)/(avg_rough*pow((avg_perem>0?avg_perem:0), 2.0/3.0));
                       beta = 5.0/3.0;
                       flux[loci][1] = locBool*alfa*pow(crossA, beta);
                       */
                       hydRadius = (avg_perem>0?crossA/avg_perem:0);
                       flux[loci][1] = locBool*sqrt(locBool*grad_y)*crossA*pow(hydRadius,2.0/3.0)/avg_rough;
                       break;
                  }
                case 2:
                  if(eletypeBool==1)
                  {
                       /* Diffusion Wave Approximation constitutive relationship: Gottardi & Venutelli, 1993 */
                       alfa = pow(pow(avg_y, 1.0/3.0),2)/(1.0*avg_rough);
                       beta = alfa;
                       flux[loci][locj] = locBool*crossA*beta*sqrt(locBool*grad_y);
                       break;
                  }
                  else
                  {
                       /*alfa = pow(pow(avg_y, 1.0/3.0),2)/(1.0*avg_rough);
                       beta = alfa;
                       flux[loci][1] = locBool*crossA*beta*sqrt(locBool*grad_y);
                       */
                       hydRadius = (avg_perem>0?crossA/avg_perem:0);
                       flux[loci][1] = locBool*sqrt(locBool*grad_y)*crossA*pow(hydRadius,2.0/3.0)/(1.0*avg_rough);
                       break;
                  }
                default:
                  if(eletypeBool==1)
                  {
                          printf("Fatal Error: Surface Overland Mode Type Is Wrong!");
                  }
                  else
                  {
                       printf("Fatal Error: River Routing Mode Type Is Wrong!");
                  }
                  exit(1);
        }
    }
}


/*    Surface Interaction between Element and River Segment    */
void OLflowFromEleToRiv(realtype sideEle_y,realtype sideEle_zmax,realtype rivX,realtype sideEleX,realtype rivY,realtype sideEleY,realtype cwr,realtype rivZmax,realtype loc_yriver,realtype **fluxriv,int loc_i,int loc_j,realtype length)
//! Computes surface flux interaction between an element and a river segment
/*! \param sideEle_y is the surface water head at the side element
    \param sideEle_zmax is surface elevation of the side element
    \param rivX is x co-ordinate of the river segment
    \param sideEleX is y co-ordinate of the side element
    \param rivY is the y co-ordinate of the river segment
    \param sideEleY is the y co-ordinate of the side element
    \param cwr is the coefficient of discharge
    \param rivZmax is the full bank elevation of the river segment
    \param loc_yriver is the river water head in the river
    \param fluxriv is the pointer to the river flux
    \param loc_i is the river segment number
    \param loc_j is 2 for left element 3 for right element
    \param length is the lenght of the river segment

*/
{
      realtype loc_cwr,dist,ele_YH,ele_Y,loc_bele;
      ele_Y = sideEle_y;
      ele_YH = sideEle_y + sideEle_zmax;
      dist = sqrt(pow(rivX - sideEleX, 2) + pow(rivY - sideEleY, 2));
      loc_cwr = cwr;
      if (rivZmax < sideEle_zmax)
      {
           loc_bele = sideEle_zmax;
      }
      else
      {
           loc_bele = rivZmax;
      }

      /*  #? where is the reference to this */
      /* ## i think its in Sorab PAnday's paper */

      if (loc_yriver > ele_YH)
      {
              if (ele_YH > loc_bele)
            {
                  fluxriv[loc_i][loc_j] = loc_cwr*2.0*sqrt(2*GRAV)*length*sqrt(loc_yriver - ele_YH)*(loc_yriver - loc_bele)/3.0;
            }
            else
            {
                if(loc_bele<loc_yriver)
                {
                      fluxriv[loc_i][loc_j] = loc_cwr*2.0*sqrt(2*GRAV)*length*sqrt(loc_yriver - loc_bele)*(loc_yriver - loc_bele)/3.0;
                     }
                else
                  {
                      fluxriv[loc_i][loc_j]=0.0;
                  }
              }

            /*fluxriv[loc_i][loc_j]=0.0; */ //uncommented last if
      }
      else{
           if (loc_yriver > loc_bele)
           {
                fluxriv[loc_i][loc_j] = -loc_cwr*2.0*sqrt(2*GRAV)*length*sqrt(ele_YH - loc_yriver)*(ele_YH - loc_bele)/3.0;
           }
           else
           {
                /* This is basicaly a lump representation */
                if(loc_bele<ele_YH)
                {
                     fluxriv[loc_i][loc_j] = -loc_cwr*2.0*sqrt(2*GRAV)*length*sqrt(ele_YH - loc_bele)*(ele_YH - loc_bele)/3.0;
                }
                else
                {
                     fluxriv[loc_i][loc_j]=0.0;
                }
           }
      }
}

/*    SubSurface Interaction between Element and River Segment    */
void GWflowFromEleToRiv(realtype sideEle_y,realtype sideEle_zmax, realtype sideEle_zmin,realtype rivX,realtype sideEleX,realtype rivY,realtype sideEleY,int loc_McPore,realtype loc_yriver,realtype loc_totyriver,realtype **fluxriv,int loc_i,int loc_j,realtype length, realtype loc_base,realtype loc_gama, realtype loc_perem,realtype loc_ksat,realtype ele_Thresh) /* delete 0.5 perimeter */
//! Computes subsurface flux interaction between an element and a river segment
/*! \param sideEle_y is the surface water head at the side element
    \param sideEle_zmax is surface elevation of the side element
    \param sideEle_zmin is bed elevation of the side element
    \param rivX is x co-ordinate of the river segment
    \param sideEleX is y co-ordinate of the side element
    \param rivY is the y co-ordinate of the river segment
    \param sideEleY is the y co-ordinate of the side element
    \param loc_McPore is the identifier for Macropore
    \param loc_yriver is the river water head
    \param loc_totyriver is the river water elevation
    \param fluxriv is the pointer to river flux
    \param loc_i is the river segment number
    \param loc_j is 4 for left element & 5 for right element
    \param length is the length of the river segment
    \param loc_base is soil base paramter : dummy
    \param loc_gama is macropore factor
    \param loc_perem is the wetted perimeter
    \param loc_ksat is the saturated hydraulic conductivity
    \param ele_Thresh is the
*/
{
     realtype loc_mpfactor,dist,ele_YH,ele_Y;

    realtype loc_sat, loc_rivK_CALIB, mp_Rzd=0.8;
    mp_Rzd = ((sideEle_zmax - sideEle_zmin - mp_Rzd) < 0)?(sideEle_zmax-sideEle_zmin):mp_Rzd;
    loc_sat = (sideEle_y-(sideEle_zmax-sideEle_zmin-mp_Rzd))/mp_Rzd;
    loc_rivK_CALIB = loc_sat>=0?(1.0+loc_sat*rivK_CALIB):1.0;

     ele_Y = sideEle_y;
     ele_YH = sideEle_y + sideEle_zmin;
     dist = sqrt(pow(rivX - sideEleX, 2) + pow(rivY - sideEleY, 2));


     if (loc_McPore == 0){
          loc_mpfactor = 1;
     }
     else if (loc_McPore == 1){
          loc_mpfactor = loc_gama;
     }

     //fluxriv[loc_i][loc_j] = length*(0.5*loc_perem)* loc_ksat *rivK_CALIB* ((loc_totyriver-ele_YH)>0? (loc_yriver>0? ((loc_totyriver-loc_yriver)>(ele_YH+ele_Thresh)?ele_Thresh:(loc_totyriver-ele_YH)):0): (ele_Y>0?((sideEle_zmin-loc_totyriver)>0?(-1.0*ele_Y):(loc_totyriver- ele_YH)):0))/dist; /* delete 0.5* perimeter */
        fluxriv[loc_i][loc_j] = length*(0.5*loc_perem)* loc_ksat *loc_rivK_CALIB* ((loc_totyriver-ele_YH)>0? (loc_yriver>0? ((loc_totyriver-loc_yriver)>(ele_YH+ele_Thresh)?ele_Thresh:(loc_totyriver-ele_YH)):0): (ele_Y>0?((sideEle_zmin-loc_totyriver)>0?(loc_totyriver- ele_YH):(loc_totyriver- ele_YH)):0))/dist;
        fluxriv[loc_i][loc_j] = fluxriv[loc_i][loc_j]>0?fluxriv[loc_i][loc_j]:fluxriv[loc_i][loc_j]*loc_mpfactor;
}
