from numpy import *


def readFile(filename):
    file_iter = open(filename, 'rU')
    dataset = []
    nameset = file_iter.readline().strip().rstrip(',').split(',')
    for line in file_iter:
        line = line.strip().rstrip(',')                         # Remove trailing comma
        temp = line.split(',')
        vec = map(float, temp)
        dataset.append(vec)
    return nameset, dataset


def randCent(dataSet, k):
    n = shape(dataSet)[1]
    centroids = mat(zeros((k, n)))  # create centroid mat
    for j in range(n):  # create random cluster centers, within bounds of each dimension
        minJ = min(dataSet[:, j])
        rangeJ = float(max(dataSet[:, j]) - minJ)
        centroids[:, j] = mat(minJ + rangeJ * random.rand(k, 1))
    return centroids


def simMeasure(vecA, vecB):
    resMat = vecA.dot(vecB.getT())
    return resMat[0,0]


def kMeans(dataSet, k, distMeas=simMeasure, createCent=randCent):
    m = shape(dataSet)[0]
    clusterAssment = mat(zeros((m, 2)))  # create mat to assign data points
    # to a centroid, also holds SE of each point
    centroids = createCent(dataSet, k)
    clusterChanged = True
    while clusterChanged:
        clusterChanged = False
        for i in range(m):  # for each data point assign it to the closest centroid
            minDist = inf
            minIndex = -1
            for j in range(k):
                distJI = distMeas(centroids[j, :], dataSet[i, :])
                if distJI < minDist:
                    minDist = distJI
                    minIndex = j
            if clusterAssment[i, 0] != minIndex:
                clusterChanged = True
            clusterAssment[i, :] = minIndex, minDist**2
        # print centroids
        for cent in range(k):  # recalculate centroids
            ptsInClust = dataSet[nonzero(clusterAssment[:, 0].A == cent)[
                0]]  # get all the point in this cluster
            # assign centroid to mean
            centroids[cent, :] = mean(ptsInClust, axis=0)
    return centroids, clusterAssment


nameset, dataset = readFile('dataset/new_case_with0.csv')
datMat = mat(dataset)
cent, cluster = kMeans(datMat, 2)
print cluster