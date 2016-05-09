# -*- coding: utf-8 -*-
"""
Created on Sat May  7 19:38:31 2016

@author: Stefan
"""
from gurobipy import *
import sys, getopt
import ast
import math

#initial interior optimizer    
def CG(n, x, populations, countries, maj1, maj2, k, pList, epsList, prec):
    
    md = Model("CG helper") 
    z = []
    ss = []
    for i in range(0,n):
        z.append(md.addVar(vtype = GRB.BINARY, name="z"+str(i)))
        
            
    z0 = md.addVar(vtype = GRB.BINARY, name="z0")
    eps = md.addVar(vtype = GRB.CONTINUOUS, name = "eps")
    md.update()    
    md.setObjective(eps, GRB.MAXIMIZE)
    for i in range(0,n):
        md.addConstr(z[i] >= 0)
    md.addConstr(eps == z0 - quicksum(z[i] * x[i] for i in range(0,n)))
    md.addConstr(quicksum(z[i] * populations[i] for i in range(0,n)) >= maj1 * z0)         
    md.addConstr(quicksum(z[i] * countries[i] for i in range(0,n)) >= maj2 * z0)

    for i in range(0,k):
        md.addConstr(eps <= epsList[i] - prec)
        md.addConstr(epsList[i] + quicksum(z[j] * pList[i][j] for j in range(0,n)) - 1 >= 0)
        
    md.update()
#    md.params.method = 2
#    md.params.crossover = 0
    md.optimize()
    
    ct = 0
    for var in md.getVars():
        if(ct < n):        
            ss.append(var.x)
        ct+=1
            
    eps = md.objVal
    for i in range(0,k):        
        assert(eps <= epsList[i])
    return (ss, eps)
    
def initialImputation(n, populations, countries, maj1, maj2):
    
    x0 = []
    for i in range(0,n):
        x0.append(1/n)
        
    epsi = -100000
    tau  = 100000
    prec = 0.000001
    
    epsList = []    
    pList = []
    
    Cr = []
    ct = 0
    
    model = Model("initialEU") 
    v = []
    
    for i in range(0,n):
        v.append(model.addVar(vtype = GRB.CONTINUOUS, name="v"+str(i)))
            
    e = model.addVar(vtype = GRB.CONTINUOUS, name="eps")
    
    model.update()    
    model.setObjective(e, GRB.MINIMIZE)
        
    for i in range(0,n):
        model.addConstr(v[i] >= 0, name="pozi"+str(i))
    
    model.addConstr(e >= 0, name="epozi")         
    model.addConstr(1 - quicksum(v[i] for i in range(0,n)) == 0, "sum1")    
    #model.params.method = 2
    #model.params.crossover = 0    
#    
    for k in range(0,n):

        if(k != 0):
            print(epsList)
            model.addConstr(e <= epsList[k-1] - prec, name="epsless"+str(k-1))
        
        model.optimize()
        status = model.status
    
        if(status != GRB.Status.OPTIMAL):
            removed = []            
            while(True):
                model.computeIIS()
                
                for c in model.getConstrs():
                    if c.IISConstr and 'pozi' not in str(c.constrName) and str(c.constrName) != 'sum1' and 'epsless' not in str(c.constrName):
                        print(c.constrName)
                        removed.append(str(c.constrName))
                        model.remove(c)
                model.optimize()
                status = model.status
                
                print(status)
                if(status == GRB.Status.OPTIMAL):
#                    ct = 0
#                    for var in model.getVars():
#                        if(ct < n):            
#                            x0[ct] = var.x
#                        ct += 1
                    break
            print("Removed constrs: ",removed)
        print("Obj value: ",model.objVal)
        for ss in Cr:
            model.addConstr(epsList[k-1] + quicksum(v[i] for i in range(0,n) if ss[i] != 0) >= \
            1 if isWinning(ss, populations, countries, maj1, maj2) else 0, "constr: "+str(ss)+" for "+str(epsList[k-1]))
        Cr = []
        
        epsi = -100000
        tau  = 100000
        model.update()
        model.optimize()
        print("Obj value now: ",model.objVal)
        while(abs(tau - epsi) > prec):   
            
            (xx, d) = CG(n, x0, populations, countries, maj1, maj2, k, pList, epsList,  prec)
            tau = d
            Sj = xx.copy()
    
            print("This is current Sj")
            print(Sj)
            if (Sj not in Cr):
                #TODO: check winning
                model.addConstr(e + quicksum(v[i] for i in range(0,n) if Sj[i] != 0) >= \
                    1 if isWinning(Sj, populations, countries, maj1, maj2) else 0, "constr"+str(Sj))
                Cr.append(Sj)
            else:
                for constr in model.getConstrs():
                    print(constr.constrName)
                print("Tau: ", tau)
                print("Eps: ", epsi)
                assert(False and "subset repeating itself")
                
            model.update()
            model.optimize()
            
            ct = 0
            for var in model.getVars():
                if(ct < n):            
                    x0[ct] = var.x
                ct += 1
            epsi = model.objVal
            
            print("Current x0, eps")
            print(x0)
            print(epsi)
            print(tau)
            print(abs(tau - epsi))
        
        for e in epsList:
            assert(epsi <= e)
            
        epsList.append(epsi)
        pList.append(x0)
        
        print("Next step!")
            
    print(pList, epsList)
    return (x0, epsi, Cr)
        
def isWinning(subset, populations, countries, maj1, maj2):
    if(sum(populations[i] for i in range(0,28) if subset[i] != 0) >= maj1 \
    and sum(countries[i] for i in range(0,28) if subset[i] !=0) >= maj2):
        return True
    return False

def main(argv):
    try:
        populations = [159, 129, 126, 119, 91, 75, 
                   39, 33, 22, 21, 20, 20, 19, 
                   19, 16, 14, 11, 10, 10, 9, 
                   8, 5, 4, 4, 3, 2, 1, 1]
        countries = [1 for i in range(0,len(populations))]
        #counts = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        
        print('No of arguments: ',len(sys.argv),' arguments')
        print('Argument list: ',str(sys.argv))
        
        try:
            opts, args = getopt.getopt(argv,"hw:q:v:z:",["weight1=","quota1=","weight2=","quota2="])
        except getopt.GetoptError:
            print('initialImputation.py -w <weight vector> -q <quota 1> -v <weight vector 2> -z <quota 2>')
            sys.exit(2)
        
        pops = []
        counts = []
        maj1 = 0
        maj2 = 0
        
        for opt, arg in opts:
            if opt == '-h':
                print('initialImputation.py -w <weight vector 1> -q <quota 1> -v <weight vector 2> -z <quota 2>')
                print('If population list is not provided, default 28 country will be used!')                
                sys.exit()
            elif opt in ("-w", "--weight1"):
                pops = arg        
            elif opt in ("-q", "--quota1"):
                maj1 = arg
            elif opt in ("-v", "--weight2"):
                counts = arg
            elif opt in ("-z", "--quota2"):
                maj2 = arg                
        
        if(pops != []):
            pops = ast.literal_eval(pops)
            maj1 = ast.literal_eval(maj1)
        else:
            pops = populations.copy()
            maj1 = 0.65 * sum(populations)
            
        if(counts != []):
            counts = ast.literal_eval(countries)
            maj2 = ast.literal_eval(maj2)
        else:
            counts = countries.copy()
            maj2 = 15
            
        n = len(pops)
        
        (x1, eps, Cs) = initialImputation(n, pops, counts, maj1, maj2)
        
        print("Obtained the following imputation:", x1)
        print("And epsilon:", eps)
        print("With added constraints:", Cs)
        
        
    except GurobiError:
            print('Error reported')
            print(GurobiError.message)
if __name__ == "__main__":
   main(sys.argv[1:])