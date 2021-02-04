def euler(x,y,f,h,x_end):
    file = open("question_one_" + str(h) + ".txt", "w+")
    while (x <= x_end):
        file.writelines([str(x)+"   ",str(y)+"\n"])
        y+=h*f(y,x)
        x+=h
    file.close()

def rk4(x,y,z,dydx,dzdx,h,x_start,x_end,flag=True):
    if (flag==True):  
        file = open("question_two_"+ str(h)+".txt", "w+")

    a=x
    b=y
    c=z

    returnlist=[[],[]]

    while (x<=x_end) or (a>=x_start):
        if (flag==True): file.writelines([str(x)+"   ",str(y)+"\n"])
        returnlist[0].append(x)
        returnlist[1].append(y)

        k1y=h*dydx(z,x)
        k1z=h*dzdx(z,y,x)
        k1b=(-h)*dydx(c,a)
        k1c=(-h)*dzdx(c,b,a)
        k2y=h*dydx(z+k1z/2,x+h/2)
        k2z=h*dzdx(z,y+k1y/2,x+h/2)
        k2b=(-h)*dydx(c+k1c/2,a+(-h)/2)
        k2c=(-h)*dzdx(c,b+k1b/2,a+(-h)/2)
        k3y=h*dydx(z+k2z/2,x+h/2)
        k3z=h*dzdx(z,y+k2y/2,x+h/2)
        k3b=(-h)*dydx(c+k2c/2,a+(-h)/2)
        k3c=(-h)*dzdx(c,b+k2b/2,a+(-h)/2)
        k4y=h*dydx(z+k3z/2,x+h)
        k4z=h*dzdx(z,y+k3y/2,x+h)
        k4b=(-h)*dydx(c+k3c/2,a+(-h))
        k4c=(-h)*dzdx(c,b+k3b/2,a+(-h))
        y+=(k1y + 2*k2y + 2*k3y + k4y)/6
        b+=(k1b + 2*k2b + 2*k3b + k4b)/6
        z+=(k1z + 2*k2z + 2*k3z + k4z)/6
        c+=(k1c + 2*k2c + 2*k3c + k4c)/6
        x+=h
        a-=h

        if (flag==True): file.writelines([str(a)+"   ",str(b)+"\n"])
    if (flag==True): file.close()
    return returnlist

def boundary(a,b,ya,yb,guess,dydx,dzdx,h,tolal):
    X=rk4(a,ya,guess,dydx,dzdx,h,a,b,False)
    Ch=guess
    i=1
    if (X[1][-1]==yb):return X
    else:
        while(X[1][-1]>=yb):
            guess-= i
            X=rk4(a,ya,guess,dydx,dzdx,h,a,b,False)
        Cl=guess
        A=rk4(a,ya,Ch,dydx,dzdx,h,a,b,False)
        B=rk4(a,ya,Cl,dydx,dzdx,h,a,b,False)
        ych=A[1][-1]
        ycl=B[1][-1]
        C=Cl+(Ch-Cl)*(yb-ycl)/(ych-ycl)
        H_L=rk4(a,ya,C,dydx,dzdx,h,a,b,False)
        yc=H_L[1][-1]
        while (abs(yc-yb)>=tolal):
            if (yc>yb) :
                Ch=C
                A=rk4(a,ya,Ch,dydx,dzdx,h,a,b,False)
                ych=A[1][-1]
                C=Cl+(Ch-Cl)*(yb-ycl)/(ych-ycl)
                H_L=rk4(a,ya,C,dydx,dzdx,h,a,b,False)
                yc=H_L[1][-1]
            else :
                Cl=C
                B=rk4(a,ya,Cl,dydx,dzdx,h,a,b,False)
                ycl=B[1][-1]
                C=Cl+(Ch-Cl)*(yb-ycl)/(ych-ycl)
                H_L=rk4(a,ya,C,dydx,dzdx,h,a,b,False)
                yc=H_L[1][-1]
    n=len(H_L[0])
    file = open("question_three_"+ str(h)+".txt", "w+")
    for i in range(n):
        file.writelines([str(H_L[0][i])+"   ",str(H_L[1][i])+"\n"])
    file.close()
    return H_L