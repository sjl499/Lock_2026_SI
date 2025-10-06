#SJL 5/17
#scripts to calculate the gradient of functions

import numpy as np

###########################################################
###########################################################
###########################################################
#define some general functions
#Function to calculate the numerical gradient of unevenly spaced points to 2nd order
#flag_calc: 0: Do not calc that differential. 1: Calc backward difference when possible. 2: Calc centered difference when possible, 3: Calc forward difference
def gradient2(x,y,flag_calc=[2,0], flag_order=[2,2]):
    
    #need to make sure the input flags are in n[ arrays for ...reasons]
    flag_calc=np.asarray(flag_calc)
    flag_order=np.asarray(flag_order)


    #check whether there are enough points for the order requested, and suitable options
    if np.all(flag_calc==0):
        print('here')
        raise Warning('No differential output requested. Returning nan')
        return np.nan*np.ones(np.size(x))
    if np.any(flag_calc)>(np.size(x)-1):
        raise ValueError('Order of differentiation requested from gradient2 function too high for number of given x-y points')
    if (flag_order[0]==1)&(flag_calc[0]==2):
        raise Warning('1st order central difference for dydx not implemented. Defaulting to forward difference')
        flag_calc[0]=1
    if (flag_calc[1]!=0)&(flag_order[1]==1):
        raise Warning('Cannot calculate 1st order d2ydx2. Defaulting to second order')
        flag_order[1]=2
    if (flag_calc[1]!=0)&(flag_calc[1]==3):
        raise Warning('Forward difference is not coded for dy2dx2. Defaulting to centered difference')
        flag_calc[1]==2
        

    
    #define the difference between all the points
    dx=np.diff(x)
    #define the output array
    dydx=np.ones(np.size(y))*np.nan

    if flag_calc[0]!=0:
        if (flag_order[0]==1)&(flag_calc[0]==1):
            dydx=gradient_backward1(x,y)
        elif (flag_order[0]==1)&(flag_calc[0]==3):
            dydx=gradient_forward1(x,y)
        elif (flag_order[0]==2):
    
            #use forward diff for first point
            dydx[0]=(y[2]*(dx[0]**2)-y[1]*((dx[1]+dx[0])**2)-y[0]*((dx[0]**2)-((dx[1]+dx[0])**2)))/\
                (dx[0]*(dx[0]+dx[1])*-dx[1])
            #use backward difference for the last point
            dydx[-1]=(y[-3]*(dx[-1]**2)-y[-2]*((dx[-1]+dx[-2])**2)-y[-1]*((dx[-1]**2)-\
                                                                          ((dx[-1]+dx[-2])**2)))/\
                (dx[-1]*(dx[-1]+dx[-2])*dx[-2])
            
            if(flag_calc[0]==1):
                #use centered difference for the second point
                dydx[1]=(y[2]*(dx[0]**2.0)-y[0]*(dx[1]**2.0)-\
                            y[1]*((dx[0]**2.0)-(dx[1]**2.0)))/\
                            (dx[0]*dx[1]*(dx[1]+dx[0]))
                    
                #use backward difference for the other points
                dydx[2:-1]=(y[0:-3]*(dx[1:-1]**2)-y[1:-2]*((dx[1:-1]+dx[0:-2])**2)-\
                            y[2:-1]*((dx[1:-1]**2)-((dx[1:-1]+dx[0:-2])**2)))/\
                        (dx[1:-1]*(dx[1:-1]+dx[0:-2])*dx[0:-2])
            elif (flag_calc[0]==2):
                #use centered difference for the middle points
                dydx[1:-1]=(y[2:]*(dx[0:-1]**2.0)-y[0:-2]*(dx[1:]**2.0)-\
                            y[1:-1]*((dx[0:-1]**2.0)-(dx[1:]**2.0)))/\
                            (dx[0:-1]*dx[1:]*(dx[1:]+dx[0:-1]))
            elif (flag_calc[0]==3):
                #use centered difference for the second to last point
                dydx[-2]=(y[-1]*(dx[-2]**2.0)-y[-3]*(dx[-1]**2.0)-\
                            y[-2]*((dx[-2]**2.0)-(dx[-1]**2.0)))/\
                            (dx[-2]*dx[-1]*(dx[-1]+dx[-2]))
                    
                #use forward diff for middle points
                dydx[0:-2]=(y[2:0]*(dx[0:-1]**2)-y[1:-1]*((dx[1:]+dx[0:-1])**2)-\
                            y[0:-2]*((dx[0:-1]**2)-((dx[1:]+dx[0:-1])**2)))/\
                    (dx[0:-1]*(dx[0:-1]+dx[1:])*-dx[1:])

            
        if flag_calc[1]==0:
            return dydx
        
    if flag_calc[1]!=0:   

            
        d2ydx2=np.ones(np.size(y))*np.nan
        A=np.ones(np.size(y))*np.nan
        B=np.ones(np.size(y))*np.nan
        C=np.ones(np.size(y))*np.nan

        #forward difference for first points
        A[0]=-2./(dx[0]*dx[1])
        B[0]=2./((dx[0]+dx[1])*dx[1])
        C[0]=2./dx[0]/(dx[0]+dx[1])
        d2ydx2[0]=A[0]*y[1]+B[0]*y[2]+C[0]*y[0]
        
        #backward difference for last point
        A[-1]=-2./(dx[-1]*dx[-2])
        B[-1]=2./((dx[-1]+dx[-2])*dx[-2])
        C[-1]=2./(dx[-1]*(dx[-1]+dx[-2]))
        d2ydx2[-1]=A[-1]*y[-2]+B[-1]*y[-3]+C[-1]*y[-1]
        
        if flag_calc[1]==1:
            
            #central difference for second point
            A[1]=-2./(-dx[0]*(dx[1]+dx[0]))
            B[1]=2./(dx[1]*(dx[1]+dx[0]))
            C[1]=2./(-dx[0]*dx[1])
            d2ydx2[1]=A[1]*y[0]+B[1]*y[2]+C[1]*y[1]
            
            #backward difference for rest of points
            A[2:]=-2./(dx[1:]*dx[0:-1])
            B[2:]=2./((dx[1:]+dx[0:-1])*dx[0:-1])
            C[2:]=2./(dx[1:]*(dx[1:]+dx[0:-1]))
            d2ydx2[2:]=A[2:]*y[1:-1]+B[2:]*y[0:-2]+C[2:]*y[2:]
            
        elif flag_calc[1]==2:
            #central difference for middle points
            A[1:-1]=-2./(-dx[0:-1]*(dx[1:]+dx[0:-1]))
            B[1:-1]=2./(dx[1:]*(dx[1:]+dx[0:-1]))
            C[1:-1]=2./(-dx[0:-1]*dx[1:])
            d2ydx2[1:-1]=A[1:-1]*y[0:-2]+B[1:-1]*y[2:]+C[1:-1]*y[1:-1]
        
        return dydx, d2ydx2


        

#Function to calculate the numerical gradient of unevenly spaced points using forward difference
def gradient_forward1(x,y):
    #define the difference between all the points
    dx=np.diff(x)
    #define the output array
    dydx=np.empty(np.size(y))

    #forward difference the first points
    dydx[0:-1]=(y[1:]-y[:-1])/dx
    
    #use backwards difference for the last point (the grad is the same as the penultimate point)
    dydx[-1]=(y[-1]-y[-2])/dx[-1]
    
    return dydx

#Function to calculate the numerical gradient of unevenly spaced points using forward difference
def gradient_backward1(x,y):
    #define the difference between all the points
    dx=np.diff(x)
    #define the output array
    dydx=np.empty(np.size(y))

    #forward difference the first point
    dydx[0]=(y[1]-y[0])/dx[0]
    
    #use backwards difference for the last point (the grad is the same as the penultimate point)
    dydx[1:]=(y[1:]-y[0:-1])/dx
    
    return dydx
