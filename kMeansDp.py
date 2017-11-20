from numpy import shape, mat, zeros, random, nonzero, mean, NINF, delete
import time
from multiprocessing import Pool


def readFile(filename):
    file_iter = open(filename, 'rU')
    dataset = []
    for line in file_iter:
        line = line.strip().rstrip(',')                         # Remove trailing comma
        temp = line.split(',')
        vec = map(float, temp)
        # print vec
        dataset.append(vec)
    print "Data loaded"
    return dataset


def getGeneName(filename):
    f = open(filename, 'r')
    genenameset = f.readline().strip().rstrip(',').split(',')
    f.close()
    return genenameset


def randCent(dataSet, k):
    n = shape(dataSet)[1]
    # print "choose representative"
    centroids = mat(zeros((k, n)))  # create centroid mat
    rprtt = 8560 * random.rand(k, 1)
    # print rprtt
    for j in range(k):  # create random cluster centers, within bounds of each dimension
        # minJ = min(dataSet[:, j])
        # rangeJ = float(max(dataSet[:, j]) - minJ)
        # centroids[:, j] = mat(minJ + rangeJ * random.rand(k, 1))
        centroids[j] = dataSet[int(rprtt[j, 0])]
    print "finish random %d cluster" % k
    return centroids


def simMeasure(vecA, vecB):
    resMat = vecA.dot(vecB.getT())
    # print resMat
    return resMat[0, 0]


def kMeans(dataSet, k, distMeas=simMeasure, createCent=randCent):
    print "enter kmeans"
    m = shape(dataSet)[0]  # line number = gene number
    clusterAssment = mat(zeros((m, 2)))  # create mat to assign data points
    centroids = createCent(dataSet, k)
    clusterChanged = True
    numIter = 0
    print "enter loop"
    while clusterChanged:
        numIter += 1
        # print numIter
        clusterChanged = False
        for i in range(m):  # assign each data point to cluster
            maxProduct = NINF
            maxIndex = -1
            for j in range(k):
                distJI = distMeas(centroids[j, :], dataSet[i, :])
                if distJI > maxProduct:
                    maxProduct = distJI
                    maxIndex = j
            if clusterAssment[i, 0] != maxIndex:
                clusterChanged = True
            clusterAssment[i, :] = maxIndex, maxProduct
        # print centroids
        for cent in range(k):  # recalculate centroids
            memberList = nonzero(clusterAssment[:, 0].A == cent)[0]
            if memberList != []:
                # get all the point in this cluster
                ptsInClust = dataSet[memberList]
                # assign centroid to mean
                centroids[cent, :] = mean(ptsInClust, axis=0)
    print numIter
    return centroids, clusterAssment


def runKMeansWithK(start, end):
    SSet = {}
    DSet = {}
    SDSet = {}
    # genenameset = getGeneName('dataset/transposed_case_0.csv')
    dataset = readFile('dataset/processed_case.csv')
    datMat = mat(dataset)
    for k in range(start, end + 1):
        print "k = %d" % k
        start_time = time.time()
        rprtts, cluAsg = kMeans(datMat, k)
        print time.time() - start_time, " seconds"
        print "cluster internal average sim"
        # calcute S
        innerAvgSim = mat(zeros((k, 2)))
        emptyClusters = []
        for i in range(k):
            memlist = nonzero(cluAsg[:, 0].A == i)[0]
            if memlist == []:       # this cluster has no member
                emptyClusters.append(i)
                continue
            print memlist
            memNum = len(memlist)
            inClusterSimSum = 0
            for m in range(memNum):
                for n in range(m + 1, memNum):
                    inClusterSimSum += simMeasure(datMat[m], datMat[n])
            numPairs = (memNum**2 - memNum) / 2
            inClusterSimSum_avg = inClusterSimSum / numPairs
            innerAvgSim[i, 0] = i
            innerAvgSim[i, 1] = inClusterSimSum_avg
        print "empty clusters", emptyClusters
        processed_innerAvgSim = delete(innerAvgSim, (emptyClusters), axis=0)
        SSet[k] = mean(processed_innerAvgSim[:, 1])
        # calculate D
        # print rprtts
        numPairs = (k**2 - k) / 2
        productSum = 0
        for i in range(k):
            for j in range(i + 1, k):
                product = simMeasure(rprtts[i], rprtts[j])
                productSum += product
        DSet[k] = productSum / numPairs
        SDSet[k] = SSet[k] / DSet[k]
    fp = open('result%d.txt' % start, 'w')
    fp.write(str(SSet))
    fp.write(str(DSet))
    fp.write(str(SDSet))
    fp.close()
    print "SSet"
    print SSet
    print "DSet"
    print DSet
    print "SDSet"
    print SDSet


p = Pool()
for x in range(33):
    p.apply_async(runKMeansWithK, args=(x * 3 + 2, x * 3 + 4,))
print 'Waiting for all process done...'
p.close()
p.join()
print 'All Done'
