import pylab as plt
import numpy as np
import os

proteins_of_interest = ['Bid_0','Bax_0','Bak_0','Bcl2_0','BclxL_0','Mcl1_0','Bad_0']
vals = np.hstack((np.linspace(0,.9,9),np.logspace(0,1,10)))
counter = 1
all = None
matrix = None
for i in proteins_of_interest:
    matrix = np.zeros((len(vals)*counter,len(vals)))
    for j in proteins_of_interest:
        filename = "fit_values_%s_%s.txt" % (i,j)
        if  os.path.exists(filename):
            tmp = np.loadtxt(filename)
            print np.shape(tmp)
            if matrix == None:
                matrix = np.loadtxt(filename)
            else:
                print np.shape(matrix)
                matrix = np.vstack((matrix,tmp))
            plot = True
            if plot == True:
                plt.imshow(tmp,interpolation='nearest',origin='lower',vmin=-1,vmax=1,cmap=plt.get_cmap('bwr'))
                plt.xticks(range(len(vals)),np.round(vals,decimals=2),rotation='vertical')
                plt.yticks(range(len(vals)),np.round(vals,decimals=2))
                plt.xlabel(j)
                plt.ylabel(i)
                plt.colorbar()
                plt.title("%s_%s" % (i,j))
                plt.savefig("fit_values_%s_%s.png" % (i,j))
                #plt.show()
                plt.clf() 
    if all == None:
        all = matrix
    else:
        all = np.hstack(((all,matrix)))
    counter +=1

all = all.T
all += all.T
print np.shape(all)
plt.imshow(all,interpolation='nearest',vmin=-1,vmax=1,origin='lower',cmap=plt.get_cmap('bwr'))
plt.xticks(np.linspace(0,len(all),len(proteins_of_interest)),proteins_of_interest,rotation='vertical')
plt.yticks(np.linspace(0,len(all),len(proteins_of_interest)),proteins_of_interest)
plt.colorbar()
plt.savefig('fit_values_all_compare.png')
plt.show()
            