/********************************************************************************
 * File        : is_sm_et.c                                                     *
 * Function    : for calculation of canopy ET, IS and snow melt       	        *
 * Version     : Nov, 2007 (2.0)                                                *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                   *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)              *
 *------------------------------------------------------------------------------*
 *                                                                              *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0.............................*
 * a) Modification of ET components from canopy, ground and transpiration       *
 * b) Addition of snow melt and throughflow processes                           *
 * c) Fully Coupled Inclusion of ET from transpiration and OVL 		        *
 * d) Rain/Snow Fraction calculation according to USACE (1956)		        *
 * e) Change in maximum interception storage due to snow accretion has been     *
 *	accounted for								*
 * f) Incorporation of Interception storage for rainfall as well as snow	*
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"      

#define EPSILON 0.05
#define UNIT_C 1440     /* 60*24 for calculation of yDot in m/min units while forcing is in m/day. */
#define multF1	1
#define multF2	1
#define multF3	1
#define multF4	1.0
realtype Interpolation(TSD *Data, realtype t);

void is_sm_et(realtype t, realtype stepsize, void *DS,N_Vector VY)
	{
  	int i;
  	realtype totEvap;
  	realtype Delta, Gamma;
  	realtype Rn, G, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a;
  	realtype isval=0,etval=0;
  	realtype fracSnow,snowRate,MeltRateGrnd,MeltRateCanopy,eltRate,MF,Ts=-3.0,Tr=1.0,To=0.0,ret;
  
  	Model_Data MD;
  
  	MD = (Model_Data)DS;
 
 	stepsize=stepsize/UNIT_C; 
  	for(i=0; i<MD->NumEle; i++)
  		{
		/* Note the dependence on physical units */
    		MD->ElePrep[i] = Interpolation(&MD->TSD_Prep[MD->Ele[i].prep-1], t);
    		Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn-1], t);
		//    G = Interpolation(&MD->TSD_G[MD->Ele[i].G-1], t);
    		T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp-1], t);
    		Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel-1], t);
    		RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity-1], t);
    		VP = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure-1], t);
    		P = 101.325*pow(10,3)*pow((293-0.0065*MD->Ele[i].zmax)/293,5.26);
    		LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC-1], t);
    		MF = multF2*Interpolation(&MD->TSD_MeltF[MD->Ele[i].meltF-1], t);	
		/******************************************************************************************/
		/*			    Snow Accumulation/Melt Calculation				  */
		/******************************************************************************************/
    		fracSnow = T<Ts?1.0:T>Tr?0:(Tr-T)/(Tr-Ts);
    		snowRate = fracSnow*MD->ElePrep[i];
		/* EleSnowGrnd, EleSnowCanopy, EleISsnowmax, MeltRateGrnd,MeltRateCanopy are the average value prorated over the whole elemental area */
    		MD->EleSnowGrnd[i]=MD->EleSnowGrnd[i]+(1-MD->Ele[i].VegFrac)*snowRate*stepsize;
    		MD->EleSnowCanopy[i]=MD->EleSnowCanopy[i]+MD->Ele[i].VegFrac*snowRate*stepsize;
		MD->EleISsnowmax[i]=MD->EleSnowCanopy[i]>0?0.003*LAI*MD->Ele[i].VegFrac:0;
		MD->EleISsnowmax[i]=multF1*MD->EleISsnowmax[i];
		if(MD->EleSnowCanopy[i]>MD->EleISsnowmax[i])
			{
			MD->EleSnowCanopy[i]=MD->EleISsnowmax[i];
			MD->EleSnowGrnd[i]=MD->EleSnowGrnd[i]+MD->EleSnowCanopy[i]-MD->EleISsnowmax[i];
			}
    		MeltRateGrnd=MeltRateCanopy=(T>To?(T-To)*MF:0);		/* Note the units for MF. */
    		if(MD->EleSnowGrnd[i]>MeltRateGrnd*stepsize)
    			{
    			MD->EleSnowGrnd[i]=MD->EleSnowGrnd[i]-MeltRateGrnd*stepsize;
    			}
    		else
    			{
    			MeltRateGrnd=MD->EleSnowGrnd[i]/stepsize;
    			MD->EleSnowGrnd[i]=0;    	
    			} 
                if(MD->EleSnowCanopy[i]>MeltRateCanopy*stepsize)
                        {
                        MD->EleSnowCanopy[i]=MD->EleSnowCanopy[i]-MeltRateCanopy*stepsize;
                        }
                else
                        {
                        MeltRateCanopy=MD->EleSnowCanopy[i]/stepsize;
                        MD->EleSnowCanopy[i]=0;
                        }
		/************************************************************************/
		/*		ThroughFall and Evaporation from canopy			*/
		/************************************************************************/
		/* EleIS, EleET[0] and ret are prorated for the whole element. Logistics are simpler if assumed in volumetric form by multiplication of Area on either side of equation*/
		MD->EleISmax[i] = multF1*MD->ISFactor[MD->Ele[i].LC-1]*LAI*MD->Ele[i].VegFrac;
		/* Note the dependence on physical units */
    		if(LAI>0.0)
    			{	 
	    		Delta = 2503*pow(10,3)*exp(17.27*T/(T+237.3))/(pow(237.3 + T, 2));
	    		Gamma = P*1.0035*0.92/(0.622*2441);
	/*    		zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
	    		cnpy_h = zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));
		    	if(LAI<2.85)	
	    			{
	    			rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
	    			}
	    		else
	    			{
	    			rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h));
	    			}
	*/	    	rl=Interpolation(&MD->TSD_RL[MD->Ele[i].LC-1], t);
	    		r_a = log(MD->Ele[i].windH/rl)*log(10*MD->Ele[i].windH/rl)/(Vel*0.16);
	    		MD->EleET[i][0] = MD->pcCal.Et0*MD->Ele[i].VegFrac*(LAI/MD->Ele[i].LAImax)*(pow((MD->EleIS[i]<0?0:MD->EleIS[i])/MD->EleISmax[i],2.0/3.0))*(Rn*(1-MD->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a))/(1000*2441000.0*(Delta+Gamma));
	    		MD->EleTF[i]=MD->EleIS[i]<=0?0:5.65*pow(10,-2)*MD->EleISmax[i]*exp(3.89*(MD->EleIS[i]<0?0:MD->EleIS[i])/MD->EleISmax[i]); /* Note the dependece on physical units*/
			MD->EleTF[i]=multF3*MD->EleTF[i];
			}
     		else
     			{
     			MD->EleET[i][0]=0.0;
			MD->EleTF[i]=0.0;
     			}
    		if(MD->EleIS[i] >= MD->EleISmax[i])
    			{
      			if(((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)>=MD->EleET[i][0]+MD->EleTF[i])
      				{
      				MD->EleETloss[i] = MD->EleET[i][0];
				ret = MD->EleTF[i]+MD->EleIS[i]-MD->EleISmax[i]+(((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)-(MD->EleET[i][0]+MD->EleTF[i]));
				MD->EleIS[i]=MD->EleISmax[i];
      				}
     			else if((((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)<MD->EleET[i][0]+MD->EleTF[i])&&(MD->EleIS[i]+stepsize*((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy-MD->EleET[i][0]-MD->EleTF[i])<=0))
     				{
     				MD->EleET[i][0]=(MD->EleET[i][0]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy));
				ret =(MD->EleTF[i]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy));  
				MD->EleIS[i] = 0;  
     				MD->EleETloss[i] =MD->EleET[i][0];
     				}
     			else
     				{
     				MD->EleIS[i] = MD->EleIS[i]+stepsize*(((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i]);
				ret =  MD->EleTF[i];   
     				MD->EleETloss[i] =MD->EleET[i][0];
  				}
    			}
    		else if((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i])*stepsize) >= MD->EleISmax[i]))
    			{
      			MD->EleETloss[i] =  MD->EleET[i][0];
      			MD->EleIS[i] = MD->EleISmax[i];
      			ret =  MD->EleTF[i]+((MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i])*stepsize)- MD->EleISmax[i]);
    			}
    		else if((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i])*stepsize) <=0))
    			{
      			if((MD->EleET[i][0]>0)||(MD->EleTF[i]>0))
				{
       	 			MD->EleET[i][0] = (MD->EleET[i][0]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy));
        			ret =(MD->EleTF[i]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)); 
				}
      			else
				{
				MD->EleET[i][0] = 0;
				ret = 0;
				}
      			MD->EleETloss[i] =  MD->EleET[i][0];
      			MD->EleIS[i] = 0;     
    			}    
    		else
    			{
      			MD->EleIS[i] = MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*MD->Ele[i].VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i])*stepsize;
      			MD->EleETloss[i] = MD->EleET[i][0];
       			ret =  MD->EleTF[i]; 
    			}
    		MD->EleNetPrep[i] = (1-MD->Ele[i].VegFrac)*(1-fracSnow)*MD->ElePrep[i]+ret+MeltRateGrnd;
		MD->EleTF[i]=ret;
//    		MD->EleNetPrep[i] = MD->ElePrep[i];
  		}
	}

