from sympy.matrices import Matrix,det,zeros
from sympy import Symbol,sin,cos,simplify,Function,solve,sqrt
from math import pi

def Christoffel1_s(g,c,a,b):
    return 0.5*(g[c,a].diff(Symbol(x[b])) +
                g[c,b].diff(Symbol(x[a])) -
                g[a,b].diff(Symbol(x[c])))
def Christoffel1(g,c,a,b):
    global christoffel1
    try:
        return christoffel1[c,a,b]
    except:
        christoffel1[c,a,b] = Christoffel1_s(g,c,a,b)
        return christoffel1[c,a,b]
def Christoffel2_s(g,I,k,l):
    global ginv
    #ginv = g.inv()
    s = 0
    n = g.shape[0]
    for m in range(n):
        s = s + ginv[I,m]*Christoffel1(g,m,k,l)
    return s
def Christoffel2(g,I,k,l):
    global christoffel2
    try:
        return christoffel2[I,k,l]
    except:
        christoffel2[I,k,l] = Christoffel2_s(g,I,k,l)
        return christoffel2[I,k,l]
def RiemannCurvature_s(g,Rho,sigma,mu,nu):
    val1 = Christoffel2(g,Rho,nu,sigma).diff(Symbol(x[mu])) - \
           Christoffel2(g,Rho,mu,sigma).diff(Symbol(x[nu]))
    val2 = 0
    val3 = 0
    n = g.shape[0]
    for lam in range(n):
        val2 = val2 + Christoffel2(g,Rho,mu,lam)*Christoffel2(
            g,lam,nu,sigma)
        val3 = val3 + Christoffel2(g,Rho,nu,lam)*Christoffel2(
            g,lam,mu,sigma)
    return (val1 + val2 - val3)
def RiemannCurvature(g,Rho,sigma,mu,nu):
    global riemanncurvature
    try:
        return riemanncurvature[Rho,sigma,mu,nu]
    except:
        riemanncurvature[Rho,sigma,mu,nu] = RiemannCurvature_s(g,Rho,sigma,mu,nu)
        return riemanncurvature[Rho,sigma,mu,nu]
        
def RicciTensor_s(g,i,j):
    s = 0
    n = g.shape[0]
    for k in range(n):
        s = s + RiemannCurvature(g,k,i,k,j)
    return s
def RicciTensor(g,i,j):
    global riccitensor
    try:
        return riccitensor[i,j]
    except:
        riccitensor[i,j] = RicciTensor_s(g,i,j)
        return riccitensor[i,j]

def ScalarCurvature(g):
    global ginv
    #ginv = g.inv()
    s = 0
    n = g.shape[0]
    for i in range(n):
        for j in range(n):
            s = s + ginv[i,j]*RicciTensor(g,i,j)
    return s

def LineElement(g, dx):
    s = 0
    n = g.shape[0]
    for mu in range(n):
        for nu in range(n):
             s = s + g[mu,nu]*Symbol(dx[mu])*Symbol(dx[nu])
    return s

def EinsteinTensor(g,i,j):
    R = ScalarCurvature(g)
    return (RicciTensor(g,i,j) - 0.5*g[i,j]*R)

def Base(i,base, bits):
    L = []
    j = 0
    n = 0
    ii = i
    for j in range(bits):
        a = ii%base
        ii = int(ii/base)
        n = n + a*base**j
        L.append(a)
    L.reverse()
    return L

def FindScalarCurvature(reset=True):
    global x,X,ginv
    global christoffel1, christoffel2, riemanncurvature, riccitensor
    if reset:
        christoffel1 = {}
        christoffel2 = {}
        riemanncurvature = {}
        riccitensor = {}
    a = 1
    if a==0:
        x = ["t","x"]
        dx = ["dt","dx"]
        n = len(x)
        [t,x1] = list(map(Symbol,x))
        Xs = ["f1","f2"]
        [f1,f2] = list(map(Function,Xs))
        X = [f1(t,x1), f2(t,x1)]
        
        X = [cos(t),x1**2*t]
    elif a==1:
        x = ["r","theta","phi"]
        dx = ["dr","dtheta","dphi"]
        n = len(x)
        [x1,x2,x3] = list(map(Symbol,x))
    
        X = [x1*sin(x2)*cos(x3),x1*sin(x2)*sin(x3),x1*cos(x2)]

    print("\nX=%s" % (str(X)))
    N = len(X)

    J = zeros( len(x), len(X))
    for i in range(N):
        for j in range(n):
            J[i,j] = X[i].diff(Symbol(x[j]))
    J = Matrix(J)
    for i in range(n**2):
        L = Base(i,n,2)
        a = L[0]
        b = L[1]
    g = (J.T)*(J)
    for i in range(N):
        for j in range(n):
            g[i,j] = g[i,j].simplify()
    s = "\nCalculating with metric tensor g_ij = \n"
    print(s)
    ginv = g.inv() 

    for i in range(n**2):
        L = Base(i,n,2)
        a = L[0]
        b = L[1]
        s = "%s%s: %s" % (str(a),str(b),str(g[a,b]))
        print(s)

    print("\n(Christoffel Symbols 1st) \n\nGamma1_ijk = ", end=' ')
    for i in range(n**3):
        L = Base(i,n,3)
        c = L[0]
        a = L[1]
        b = L[2]
        print(str(c)+str(a)+str(b)+":",Christoffel1(g,c,a,b),",", end=' ')

    print("\n\n(Christoffel Symbols 2nd) Gamma2_Ijk = ", end=' ')
    for i in range(n**3):
        L = Base(i,n,3)
        c = L[0]
        a = L[1]
        b = L[2]
        print(str(c)+str(a)+str(b)+":",Christoffel2(g,c,a,b),",", end=' ')


    
    dA = sqrt(det(g))
    for j in range(n):
        dA *= Symbol(dx[j])
    s = "\ndA = sqrt(det(g[i,j]))*dx*dy = %s" % (str(dA))
    print(s)
    
    s = "\n\nds**2 = %s" % str(LineElement(g,dx))
    print(s)


    s = "\n\n(Ricci Tensor) R_ij = \n"
    print(s)
    for i in range(n**2):
        L = Base(i,n,2)
        a = L[0]
        b = L[1]
        R2 = RicciTensor(g,a,b).expand().simplify()
        s = "%s%s: %s" % (str(a),str(b),str(R2))
        print(s)
        
    R = ScalarCurvature(g).expand().simplify()
    s = "\n\nScalar Curvature R = %s" % str(R)
    print(s)
    return

FindScalarCurvature()

