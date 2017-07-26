import numpy as np
pyseed = 1
np.random.seed(pyseed)

## ********************************************************** ##
## ********************  Some functions  ******************** ##
## ********************************************************** ##
#  to be sorted
#  to read, skip to the end of this session and continue

## Create a matrix to determine the parameters (e.g. conductance) of the cell-mechanism
#  For N-cells and M-parameters, define the matrix to be M-by-N matrix
#
#  On the way, loop through the cell and assign to them.

HetIndexHermann2007 = {0:'rho_NaV', \
                                        1:'rho_NaK', \
                                        2:'rho_NCX', \
                                        3:'rho_KATP', \
                                        4:'rho_KV', \
                                        5:'rho_sKCa', \
                                        6:'rho_KCa', \
                                        7:'rho_CaL', \
                                        8:'rho_CaT', \
                                        9:'rho_PMCA'}

HetDictHermann2007 = {'rho_NaV':(1.4, 0.05*1.4), \
                                        'rho_NaK':(2000., 0.05*2000.), \
                                        'rho_NCX':(14., 0.05*14.), \
                                        'rho_KATP':(0.08, 0.05*0.08), \
                                        'rho_KV':(4.9, 0.05*4.9), \
                                        'rho_sKCa':(0.6, 0.05*0.6), \
                                        'rho_KCa':(0.4, 0.05*0.4), \
                                        'rho_CaL':(0.7, 0.05*0.7), \
                                        'rho_CaT':(0.15, 0.05*0.15), \
                                        'rho_PMCA':(1100., 0.05*1100.)}


def setHeteroHermann2007(b_cell,HetMatrix,i):
    # Passing b_cell, HetMatrix as a mutable object
    # b_cell: neuron h object constants function soma(0.5)
    # HetMatrix: matrix to store heterogeneity of cells
    # i: index (id) of b_cell
    b_cell.soma(0.5).rho_NaV_jna = HetMatrix[0,i] = np.random.normal(HetDictHermann2007['rho_NaV'][0],HetDictHermann2007['rho_NaV'][1])
    b_cell.soma(0.5).gbar_nap = HetMatrix[0,i]*b_cell.g_NaV_bar

    b_cell.soma(0.5).rho_NaK_jna = HetMatrix[1,i] = np.random.normal(HetDictHermann2007['rho_NaK'][0],HetDictHermann2007['rho_NaK'][1])
    b_cell.soma(0.5).rho_NaK_jk = HetMatrix[1,i]
    b_cell.soma(0.5).gbar_nak = HetMatrix[1,i]*b_cell.I_NaK_bar

    b_cell.soma(0.5).rho_NCX_jna = HetMatrix[2,i] = np.random.normal(HetDictHermann2007['rho_NCX'][0],HetDictHermann2007['rho_NCX'][1])
    b_cell.soma(0.5).rho_NCX_jca = HetMatrix[2,i]
    b_cell.soma(0.5).gbar_ncx = HetMatrix[2,i]*b_cell.I_NCX_bar

    b_cell.soma(0.5).rho_KATP_jk = HetMatrix[3,i] = np.random.normal(HetDictHermann2007['rho_KATP'][0],HetDictHermann2007['rho_KATP'][1])
    b_cell.soma(0.5).gbar_katp = HetMatrix[3,i]*b_cell.g_KATP_bar

    b_cell.soma(0.5).rho_KV_jk = HetMatrix[4,i] = np.random.normal(HetDictHermann2007['rho_KV'][0],HetDictHermann2007['rho_KV'][1])
    b_cell.soma(0.5).gbar_kv = HetMatrix[4,i]*b_cell.g_KV_bar

    b_cell.soma(0.5).rho_sKCa_jk = HetMatrix[5,i] = np.random.normal(HetDictHermann2007['rho_sKCa'][0],HetDictHermann2007['rho_sKCa'][1])
    b_cell.soma(0.5).gbar_skca = HetMatrix[5,i]*b_cell.g_sKCa_bar

    b_cell.soma(0.5).rho_KCa_jk = HetMatrix[6,i] = np.random.normal(HetDictHermann2007['rho_KCa'][0],HetDictHermann2007['rho_KCa'][1])
    b_cell.soma(0.5).gbar_kca = HetMatrix[6,i]*b_cell.g_KCa_bar

    b_cell.soma(0.5).rho_CaL_jca = HetMatrix[7,i] = np.random.normal(HetDictHermann2007['rho_CaL'][0],HetDictHermann2007['rho_CaL'][1])
    b_cell.soma(0.5).gbar_cal = HetMatrix[7,i]*b_cell.g_CaL_bar

    b_cell.soma(0.5).rho_CaT_jca = HetMatrix[8,i] = np.random.normal(HetDictHermann2007['rho_CaT'][0],HetDictHermann2007['rho_CaT'][1])
    b_cell.soma(0.5).gbar_cat = HetMatrix[8,i]*b_cell.g_CaT_bar

    b_cell.soma(0.5).rho_PMCA_jca = HetMatrix[9,i] = np.random.normal(HetDictHermann2007['rho_PMCA'][0],HetDictHermann2007['rho_PMCA'][1])
    b_cell.soma(0.5).gbar_pmca = HetMatrix[9,i]*b_cell.I_PMCA_bar

    b_cell.soma(0.5).nkatp_katp = HetMatrix[10,i] = np.random.normal(-6.5,0.3)
    b_cell.soma(0.5).kappa_KATP_jk = HetMatrix[10,i]


def loadHeteroHermann2007(b_cell,HetMatrix,i):
    # Passing b_cell, HetMatrix as a mutable object
    # b_cell: neuron h object constants function soma(0.5)
    # HetMatrix: matrix to assign heterogeneity of cells
    # i: index (id) of b_cell
    b_cell.soma(0.5).rho_NaV_jna = HetMatrix[0,i]
    b_cell.soma(0.5).gbar_nap = HetMatrix[0,i]*b_cell.g_NaV_bar

    b_cell.soma(0.5).rho_NaK_jna = HetMatrix[1,i]
    b_cell.soma(0.5).rho_NaK_jk = HetMatrix[1,i]
    b_cell.soma(0.5).gbar_nak = HetMatrix[1,i]*b_cell.I_NaK_bar

    b_cell.soma(0.5).rho_NCX_jna = HetMatrix[2,i]
    b_cell.soma(0.5).rho_NCX_jca = HetMatrix[2,i]
    b_cell.soma(0.5).gbar_ncx = HetMatrix[2,i]*b_cell.I_NCX_bar

    b_cell.soma(0.5).rho_KATP_jk = HetMatrix[3,i]
    b_cell.soma(0.5).gbar_katp = HetMatrix[3,i]*b_cell.g_KATP_bar

    b_cell.soma(0.5).rho_KV_jk = HetMatrix[4,i]
    b_cell.soma(0.5).gbar_kv = HetMatrix[4,i]*b_cell.g_KV_bar

    b_cell.soma(0.5).rho_sKCa_jk = HetMatrix[5,i]
    b_cell.soma(0.5).gbar_skca = HetMatrix[5,i]*b_cell.g_sKCa_bar

    b_cell.soma(0.5).rho_KCa_jk = HetMatrix[6,i]
    b_cell.soma(0.5).gbar_kca = HetMatrix[6,i]*b_cell.g_KCa_bar

    b_cell.soma(0.5).rho_CaL_jca = HetMatrix[7,i]
    b_cell.soma(0.5).gbar_cal = HetMatrix[7,i]*b_cell.g_CaL_bar

    b_cell.soma(0.5).rho_CaT_jca = HetMatrix[8,i]
    b_cell.soma(0.5).gbar_cat = HetMatrix[8,i]*b_cell.g_CaT_bar

    b_cell.soma(0.5).rho_PMCA_jca = HetMatrix[9,i]
    b_cell.soma(0.5).gbar_pmca = HetMatrix[9,i]*b_cell.I_PMCA_bar

    b_cell.soma(0.5).nkatp_katp = HetMatrix[10,i]
    b_cell.soma(0.5).kappa_KATP_jk = HetMatrix[10,i]


## Define the CoupledMatrix if provided a coordinate data
#  Procedure:
#  - load the coordinate
#  - define distance matrix
#    - by loop through the coordinates twice
#  - set cutoff distance cutoffConnect ~ 20 um
#  - convert it into binary matrix

def genCoupleMatrix(coorData, cutoff=20):
    # Generate a binary CoupleMatrix using islet coordinates data
    #
    # if given 3 columns coordinate data of one cell type (e.g. b-cell)
    if isinstance(coorData, (str, unicode)):
        coorData = np.loadtxt(coorData)
    # else assume it is a numpy array type with equivalent format
    # expecting format of N-by-3 array (for N cells)
    N = np.shape(coorData)[0]        # change here if different format
    # create empty CoupleMatrix
    CoupleMatrix = np.zeros([N,N])
    for i in xrange(N):
        for j in xrange(N):
            diff = coorData[i,:] - coorData[j,:]
            CoupleMatrix[i,j] = np.inner(diff,diff)
    cutoff = cutoff**2            # compare in distance-square
    # convert it to binary matrix
    CoupleMatrix[CoupleMatrix<cutoff] = 1
    CoupleMatrix[CoupleMatrix>cutoff] = 0
    return CoupleMatrix

## Tested
#A = np.array([[1, 2, 3],[3,2,1],[10,5,1],[100,1,2]])
#A = "testCoorData.txt"
#B = genCoupleMatrix(A)


## Define array to set b-cell hugs using cell index (id)
## Define b-cell hugs properties (potentially another .hoc)


## output CoupledMatrix, parameter matrix, hugs array
def setSimOutput(b_cells):
    # Return 3 lists: time, Vm, Ca_i
    # Except time, other outputs are list of variable of all cells
    # b_cells: a list containing b_cell objects
    n = len(b_cells)
    time = h.Vector()
    time.record(h._ref_t)
    vrec = []
    carec = []
    for i in range(n):
        vrec.append(h.Vector())
        carec.append(h.Vector())
        vrec[i].record(cell[i].soma(0.5)._ref_v)
        carec[i].record(cell[i].soma(0.5)._ref_cai)
    return time, vrec, carec

def convertSimOutput(listOutput,sparse=1):
    # Convert a list of h Vector object output to numpy array object
    # DO NOT USE if just want to convert a few of the cells in the whole list
    n = len(listOutput)
    for i in range(n):
        listOutput[i] = np.array(listOutput[i])[::sparse]
    return listOutput

## ********************************************************** ##
## ********************************************************** ##

